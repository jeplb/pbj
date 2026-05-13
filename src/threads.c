#include "threads.h"
#include "merge.h"
#include "adapter.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

typedef struct {
    pbj_workspace_t   ws;
    const pbj_params_t *params;
    pbj_read_t        *r1_buf;
    pbj_read_t        *r2_buf;
    pbj_merge_result_t *out_buf;
    int                start;
    int                end;
    pthread_t          tid;
} worker_t;

static void *worker_main(void *arg) {
    worker_t *w = (worker_t*)arg;
    for (int i = w->start; i < w->end; i++) {
        pbj_try_merge(&w->r1_buf[i], &w->r2_buf[i], &w->ws, w->params, &w->out_buf[i]);
    }
    return NULL;
}

static int account_and_emit(int i,
                             pbj_read_t *r1_buf,
                             pbj_read_t *r2_buf,
                             pbj_merge_result_t *out_buf,
                             pbj_fastq_writer_t *aw,
                             pbj_fastq_writer_t *fw,
                             pbj_fastq_writer_t *rw,
                             pbj_stats_t *stats,
                             const pbj_params_t *p) {
    stats->total_pairs++;
    pbj_merge_result_t *res = &out_buf[i];
    int rc = 0;
    if (res->status == PBJ_MERGED_OK) {
        stats->merged++;
        /* merged read qualities are already phred+33 (written by merge.c). */
        if (pbj_fastq_write(aw, &res->merged, 33) != 0) rc = -1;
    } else {
        switch (res->status) {
            case PBJ_NO_CANDIDATE:    stats->no_candidate++;   break;
            case PBJ_FAIL_PVALUE:     stats->fail_pvalue++;    break;
            case PBJ_FAIL_LENGTH:     stats->fail_length++;    break;
            case PBJ_FAIL_AMBIGUOUS:  stats->fail_ambiguous++; break;
            case PBJ_FAIL_MIN_QUAL:   stats->fail_min_qual++;  break;
            case PBJ_FAIL_MAX_EE:     stats->fail_max_ee++;    break;
            case PBJ_FAIL_MONO:       stats->fail_mono++;      break;
            default:
                /* an unknown pbj_fail_* would silently break the stats-sum
                   invariant (total_pairs incremented but no sub-counter).
                   abort instead. */
                fprintf(stderr, "pbj: internal error, unknown status %d\n",
                        (int)res->status);
                abort();
        }
        if (pbj_fastq_write(fw, &r1_buf[i], p->quality_offset) != 0) rc = -1;
        if (pbj_fastq_write(rw, &r2_buf[i], p->quality_offset) != 0) rc = -1;
    }
    if (p->verbose) {
        const char *tag;
        switch (res->status) {
            case PBJ_MERGED_OK:       tag = "OK";       break;
            case PBJ_NO_CANDIDATE:    tag = "NOSEED";   break;
            case PBJ_FAIL_PVALUE:     tag = "PVAL";     break;
            case PBJ_FAIL_LENGTH:     tag = "LEN";      break;
            case PBJ_FAIL_AMBIGUOUS:  tag = "AMBIG";    break;
            case PBJ_FAIL_MIN_QUAL:   tag = "MINQUAL";  break;
            case PBJ_FAIL_MAX_EE:     tag = "EE";       break;
            case PBJ_FAIL_MONO:       tag = "MONO";     break;
            default:                  tag = "UNK";      break;
        }
        fprintf(stderr, "  %s s=%d olen=%d m=%d mm=%d nc=%d p=%.3g name=%s\n",
                tag, res->offset, res->overlap_len, res->matches, res->mismatches,
                res->n_candidates, res->p_value,
                r1_buf[i].name ? r1_buf[i].name : "");
    }
    return rc;
}

int pbj_run_pipeline(const pbj_params_t *p,
                     pbj_fastq_reader_t *r1_reader,
                     pbj_fastq_reader_t *r2_reader,
                     pbj_fastq_writer_t *aw,
                     pbj_fastq_writer_t *fw,
                     pbj_fastq_writer_t *rw,
                     pbj_stats_t *stats) {
    int N = p->threads < 1 ? 1 : p->threads;
    int B = p->batch_size < 1 ? 1 : p->batch_size;
    if (B < N) B = N;

    pbj_read_t         *r1_buf  = (pbj_read_t*)calloc((size_t)B, sizeof(pbj_read_t));
    pbj_read_t         *r2_buf  = (pbj_read_t*)calloc((size_t)B, sizeof(pbj_read_t));
    pbj_merge_result_t *out_buf = (pbj_merge_result_t*)calloc((size_t)B, sizeof(pbj_merge_result_t));
    worker_t           *workers = NULL;

    int rc = 0;
    if (!r1_buf || !r2_buf || !out_buf) {
        fprintf(stderr, "pbj: out of memory allocating batch buffers\n");
        rc = 3;
        goto done;
    }
    for (int i = 0; i < B; i++) {
        pbj_read_init(&r1_buf[i]);
        pbj_read_init(&r2_buf[i]);
        pbj_read_init(&out_buf[i].merged);
    }

    workers = (worker_t*)calloc((size_t)N, sizeof(worker_t));
    if (!workers) {
        fprintf(stderr, "pbj: out of memory allocating workers\n");
        rc = 3;
        goto done;
    }
    for (int i = 0; i < N; i++) {
        workers[i].params  = p;
        workers[i].r1_buf  = r1_buf;
        workers[i].r2_buf  = r2_buf;
        workers[i].out_buf = out_buf;
        pbj_workspace_init(&workers[i].ws, p);
    }

    int first_pair_checked = 0;
    int interleaved        = (r1_reader == r2_reader);
    for (;;) {
        int n = 0;
        for (; n < B; n++) {
            int rc1 = pbj_fastq_read(r1_reader, &r1_buf[n]);
            if (rc1 == 0) break;
            if (rc1 < 0) { fprintf(stderr, "pbj: failed to read R1 record\n"); rc = 2; goto done; }
            int rc2 = pbj_fastq_read(r2_reader, &r2_buf[n]);
            if (rc2 == 0) { fprintf(stderr, "pbj: R2 ended before R1 (input has an odd record count)\n"); rc = 2; goto done; }
            if (rc2 < 0) { fprintf(stderr, "pbj: failed to read R2 record\n"); rc = 2; goto done; }
            if (!pbj_pair_names_match(&r1_buf[n], &r2_buf[n])) {
                fprintf(stderr, "pbj: read name mismatch at pair %llu (R1='%.*s' R2='%.*s')\n",
                        (unsigned long long)(stats->total_pairs + (uint64_t)n),
                        (int)r1_buf[n].name_len, r1_buf[n].name,
                        (int)r2_buf[n].name_len, r2_buf[n].name);
                rc = 2;
                goto done;
            }
            /* on the first pair, verify direction via casava 1.8+ tokens
               when present. catches both interleaved-r2-first input and a
               swapped -f/-r pair on separate files. non-casava headers
               cannot be checked and the file order is trusted. */
            if (!first_pair_checked) {
                int d1 = pbj_casava_direction(&r1_buf[n]);
                int d2 = pbj_casava_direction(&r2_buf[n]);
                if (d1 && d2 && !(d1 == 1 && d2 == 2)) {
                    const char *what = interleaved ? "interleaved input" : "input pair";
                    fprintf(stderr,
                        "pbj: %s appears to be in (R2,R1) order "
                        "(direction tokens R1=%d, R2=%d)\n", what, d1, d2);
                    rc = 2;
                    goto done;
                }
                first_pair_checked = 1;
            }
            if (p->adapter_forward) {
                if (pbj_adapter_trim(&r1_buf[n], p->adapter_forward,
                                      p->adapter_min_match, p->adapter_match_frac)) {
                    stats->adapter_trimmed_r1++;
                }
            }
            if (p->adapter_reverse) {
                if (pbj_adapter_trim(&r2_buf[n], p->adapter_reverse,
                                      p->adapter_min_match, p->adapter_match_frac)) {
                    stats->adapter_trimmed_r2++;
                }
            }
        }
        if (n == 0) break;

        if (N == 1) {
            workers[0].start = 0;
            workers[0].end   = n;
            (void)worker_main(&workers[0]);
        } else {
            int chunk = (n + N - 1) / N;
            for (int i = 0; i < N; i++) {
                int s = i * chunk;
                int e = s + chunk;
                if (s > n) s = n;
                if (e > n) e = n;
                workers[i].start = s;
                workers[i].end   = e;
            }
            for (int i = 0; i < N; i++) {
                if (workers[i].start >= workers[i].end) continue;
                if (pthread_create(&workers[i].tid, NULL, worker_main, &workers[i]) != 0) {
                    fprintf(stderr, "pbj: pthread_create failed\n");
                    rc = 3;
                    for (int j = 0; j < i; j++) {
                        if (workers[j].start < workers[j].end) pthread_join(workers[j].tid, NULL);
                    }
                    goto done;
                }
            }
            for (int i = 0; i < N; i++) {
                if (workers[i].start < workers[i].end) pthread_join(workers[i].tid, NULL);
            }
        }

        for (int i = 0; i < n; i++) {
            if (account_and_emit(i, r1_buf, r2_buf, out_buf, aw, fw, rw, stats, p) != 0) {
                fprintf(stderr, "pbj: write error\n");
                rc = 2;
                goto done;
            }
        }
    }
    /* in non-interleaved mode, check whether r2 has unread records past r1's
       eof. in interleaved mode this is moot because the two readers are the
       same stream. */
    if (rc == 0 && !interleaved) {
        pbj_read_t tmp;
        pbj_read_init(&tmp);
        int rc2 = pbj_fastq_read(r2_reader, &tmp);
        if (rc2 == 1) {
            fprintf(stderr, "pbj: R1 ended before R2\n");
            rc = 2;
        }
        pbj_read_free(&tmp);
    }

done:
    if (workers) {
        for (int i = 0; i < N; i++) pbj_workspace_free(&workers[i].ws);
        free(workers);
    }
    if (r1_buf) {
        for (int i = 0; i < B; i++) pbj_read_free(&r1_buf[i]);
        free(r1_buf);
    }
    if (r2_buf) {
        for (int i = 0; i < B; i++) pbj_read_free(&r2_buf[i]);
        free(r2_buf);
    }
    if (out_buf) {
        for (int i = 0; i < B; i++) pbj_read_free(&out_buf[i].merged);
        free(out_buf);
    }
    return rc;
}
