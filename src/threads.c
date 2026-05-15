#include "threads.h"
#include "merge.h"
#include "adapter.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

/* persistent thread pool. each worker spins on a condvar waiting for a
   chunk of pairs; the main thread dispatches a [start,end) range per
   batch and waits on a "done" condvar. avoids pthread_create / join per
   batch (10-50us each) which becomes the dominant overhead at small
   batch sizes. */
typedef struct {
    pthread_t          tid;
    pthread_mutex_t    mu;
    pthread_cond_t     work_cv;   /* main signals worker: new chunk ready */
    pthread_cond_t     done_cv;   /* worker signals main: chunk finished  */
    pbj_workspace_t    ws;
    const pbj_params_t *params;
    pbj_read_t        *r1_buf;
    pbj_read_t        *r2_buf;
    pbj_merge_result_t *out_buf;
    int                start;
    int                end;
    int                has_work;  /* 1 while a chunk is pending or running */
    int                exit_flag; /* tells the worker to return */
    uint64_t           local_trim_r1; /* per-batch local counters, summed  */
    uint64_t           local_trim_r2; /* into stats after the batch joins  */
} worker_t;

static inline void worker_process_one(worker_t *w, int i) {
    if (w->params->adapter_forward) {
        if (pbj_adapter_trim(&w->r1_buf[i], w->params->adapter_forward,
                              w->params->adapter_min_match,
                              w->params->adapter_match_frac)) {
            w->local_trim_r1++;
        }
    }
    if (w->params->adapter_reverse) {
        if (pbj_adapter_trim(&w->r2_buf[i], w->params->adapter_reverse,
                              w->params->adapter_min_match,
                              w->params->adapter_match_frac)) {
            w->local_trim_r2++;
        }
    }
    pbj_try_merge(&w->r1_buf[i], &w->r2_buf[i], &w->ws, w->params, &w->out_buf[i]);
}

static void *worker_main(void *arg) {
    worker_t *w = (worker_t*)arg;
    for (;;) {
        pthread_mutex_lock(&w->mu);
        while (!w->has_work && !w->exit_flag) {
            pthread_cond_wait(&w->work_cv, &w->mu);
        }
        if (w->exit_flag) {
            pthread_mutex_unlock(&w->mu);
            return NULL;
        }
        int s = w->start;
        int e = w->end;
        pthread_mutex_unlock(&w->mu);

        for (int i = s; i < e; i++) {
            worker_process_one(w, i);
        }

        pthread_mutex_lock(&w->mu);
        w->has_work = 0;
        pthread_cond_signal(&w->done_cv);
        pthread_mutex_unlock(&w->mu);
    }
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

/* a batch of B paired records and their merge results. the pipeline keeps
   a small ring of these so the reader thread can fill slot N+1 while
   workers + writer are still busy with slot N. */
typedef struct {
    pbj_read_t        *r1;
    pbj_read_t        *r2;
    pbj_merge_result_t *out;
    int                n;     /* records actually read; may be < B at EOF */
} batch_t;

#define NSLOTS 2

typedef struct {
    batch_t             slots[NSLOTS];
    int                 B;
    int                 head;          /* next slot for reader to fill   */
    int                 tail;          /* next slot for main to consume  */
    int                 count;         /* number of full slots in [tail,head) */
    pthread_mutex_t     mu;
    pthread_cond_t      not_full_cv;   /* signaled when main empties a slot */
    pthread_cond_t      not_empty_cv;  /* signaled when reader fills a slot */
    int                 reader_done;   /* reader saw EOF or fatal error  */
    int                 reader_rc;     /* 0 ok, 2 i/o error              */
    int                 fatal;         /* main asks reader to stop early */
    pbj_fastq_reader_t *r1_reader;
    pbj_fastq_reader_t *r2_reader;
    int                 interleaved;
} input_ring_t;

static void *reader_main(void *arg) {
    input_ring_t *pl = (input_ring_t*)arg;
    int first_pair_checked = 0;
    uint64_t pair_counter = 0;

    for (;;) {
        pthread_mutex_lock(&pl->mu);
        while (pl->count == NSLOTS && !pl->fatal) {
            pthread_cond_wait(&pl->not_full_cv, &pl->mu);
        }
        if (pl->fatal) { pthread_mutex_unlock(&pl->mu); return NULL; }
        int idx = pl->head;
        pthread_mutex_unlock(&pl->mu);

        batch_t *b = &pl->slots[idx];
        int n = 0;
        int err = 0;
        int eof = 0;
        for (; n < pl->B; n++) {
            int rc1 = pbj_fastq_read(pl->r1_reader, &b->r1[n]);
            if (rc1 == 0) { eof = 1; break; }
            if (rc1 < 0) {
                fprintf(stderr, "pbj: failed to read R1 record\n");
                err = 2; break;
            }
            int rc2 = pbj_fastq_read(pl->r2_reader, &b->r2[n]);
            if (rc2 == 0) {
                fprintf(stderr, "pbj: R2 ended before R1 (input has an odd record count)\n");
                err = 2; break;
            }
            if (rc2 < 0) {
                fprintf(stderr, "pbj: failed to read R2 record\n");
                err = 2; break;
            }
            if (!pbj_pair_names_match(&b->r1[n], &b->r2[n])) {
                fprintf(stderr, "pbj: read name mismatch at pair %llu (R1='%.*s' R2='%.*s')\n",
                        (unsigned long long)pair_counter,
                        (int)b->r1[n].name_len, b->r1[n].name,
                        (int)b->r2[n].name_len, b->r2[n].name);
                err = 2; break;
            }
            /* on the first pair, verify direction via casava 1.8+ tokens
               when present. catches both interleaved-r2-first input and a
               swapped -f/-r pair on separate files. non-casava headers
               cannot be checked and the file order is trusted. */
            if (!first_pair_checked) {
                int d1 = pbj_casava_direction(&b->r1[n]);
                int d2 = pbj_casava_direction(&b->r2[n]);
                if (d1 && d2 && !(d1 == 1 && d2 == 2)) {
                    const char *what = pl->interleaved ? "interleaved input" : "input pair";
                    fprintf(stderr,
                        "pbj: %s appears to be in (R2,R1) order "
                        "(direction tokens R1=%d, R2=%d)\n", what, d1, d2);
                    err = 2; break;
                }
                first_pair_checked = 1;
            }
            pair_counter++;
        }
        /* on R1 EOF (eof = 1), verify R2 also reached EOF unless interleaved
           (same stream). picks up "R1 ended before R2" cases. */
        if (eof && !pl->interleaved && err == 0) {
            pbj_read_t tmp;
            pbj_read_init(&tmp);
            int rc2 = pbj_fastq_read(pl->r2_reader, &tmp);
            if (rc2 == 1) {
                fprintf(stderr, "pbj: R1 ended before R2\n");
                err = 2;
            }
            pbj_read_free(&tmp);
        }

        b->n = n;

        pthread_mutex_lock(&pl->mu);
        if (err) {
            pl->reader_rc = err;
            pl->reader_done = 1;
            pthread_cond_broadcast(&pl->not_empty_cv);
            pthread_mutex_unlock(&pl->mu);
            return NULL;
        }
        if (n > 0) {
            pl->head = (pl->head + 1) % NSLOTS;
            pl->count++;
            pthread_cond_signal(&pl->not_empty_cv);
        }
        if (eof) {
            pl->reader_done = 1;
            pthread_cond_broadcast(&pl->not_empty_cv);
            pthread_mutex_unlock(&pl->mu);
            return NULL;
        }
        pthread_mutex_unlock(&pl->mu);
    }
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

    worker_t      *workers = NULL;
    int            workers_started = 0;
    input_ring_t   ring;
    int            ring_inited = 0;
    int            reader_started = 0;
    pthread_t      reader_tid;

    int rc = 0;

    memset(&ring, 0, sizeof(ring));
    ring.B = B;
    ring.r1_reader = r1_reader;
    ring.r2_reader = r2_reader;
    ring.interleaved = (r1_reader == r2_reader);
    for (int s = 0; s < NSLOTS; s++) {
        ring.slots[s].r1  = (pbj_read_t*)calloc((size_t)B, sizeof(pbj_read_t));
        ring.slots[s].r2  = (pbj_read_t*)calloc((size_t)B, sizeof(pbj_read_t));
        ring.slots[s].out = (pbj_merge_result_t*)calloc((size_t)B, sizeof(pbj_merge_result_t));
        if (!ring.slots[s].r1 || !ring.slots[s].r2 || !ring.slots[s].out) {
            fprintf(stderr, "pbj: out of memory allocating batch buffers\n");
            rc = 3; goto done;
        }
        for (int i = 0; i < B; i++) {
            pbj_read_init(&ring.slots[s].r1[i]);
            pbj_read_init(&ring.slots[s].r2[i]);
            pbj_read_init(&ring.slots[s].out[i].merged);
        }
    }
    if (pthread_mutex_init(&ring.mu, NULL) != 0 ||
        pthread_cond_init(&ring.not_full_cv, NULL) != 0 ||
        pthread_cond_init(&ring.not_empty_cv, NULL) != 0) {
        fprintf(stderr, "pbj: pthread sync init failed\n");
        rc = 3; goto done;
    }
    ring_inited = 1;

    workers = (worker_t*)calloc((size_t)N, sizeof(worker_t));
    if (!workers) {
        fprintf(stderr, "pbj: out of memory allocating workers\n");
        rc = 3;
        goto done;
    }
    for (int i = 0; i < N; i++) {
        workers[i].params  = p;
        pbj_workspace_init(&workers[i].ws, p);
    }
    /* skip pthread setup for N==1: the main thread does the work itself.
       avoids ~20% slowdown vs an inline call when nothing is parallel. */
    if (N > 1) {
        for (int i = 0; i < N; i++) {
            if (pthread_mutex_init(&workers[i].mu, NULL) != 0 ||
                pthread_cond_init(&workers[i].work_cv, NULL) != 0 ||
                pthread_cond_init(&workers[i].done_cv, NULL) != 0) {
                fprintf(stderr, "pbj: pthread sync init failed\n");
                rc = 3; goto done;
            }
            if (pthread_create(&workers[i].tid, NULL, worker_main, &workers[i]) != 0) {
                fprintf(stderr, "pbj: pthread_create failed\n");
                rc = 3; goto done;
            }
            workers_started++;
        }
    }

    if (pthread_create(&reader_tid, NULL, reader_main, &ring) != 0) {
        fprintf(stderr, "pbj: pthread_create failed for reader\n");
        rc = 3; goto done;
    }
    reader_started = 1;

    for (;;) {
        pthread_mutex_lock(&ring.mu);
        while (ring.count == 0 && !ring.reader_done) {
            pthread_cond_wait(&ring.not_empty_cv, &ring.mu);
        }
        if (ring.count == 0) {
            if (ring.reader_rc != 0) rc = ring.reader_rc;
            pthread_mutex_unlock(&ring.mu);
            break;
        }
        int idx = ring.tail;
        pthread_mutex_unlock(&ring.mu);

        batch_t *b = &ring.slots[idx];
        int n = b->n;

        for (int i = 0; i < N; i++) {
            workers[i].r1_buf       = b->r1;
            workers[i].r2_buf       = b->r2;
            workers[i].out_buf      = b->out;
            workers[i].local_trim_r1 = 0;
            workers[i].local_trim_r2 = 0;
        }
        if (N == 1) {
            workers[0].start = 0; workers[0].end = n;
            for (int i = 0; i < n; i++) worker_process_one(&workers[0], i);
        } else {
            int chunk = (n + N - 1) / N;
            for (int i = 0; i < N; i++) {
                int s = i * chunk;
                int e = s + chunk;
                if (s > n) s = n;
                if (e > n) e = n;
                if (s >= e) continue;
                pthread_mutex_lock(&workers[i].mu);
                workers[i].start    = s;
                workers[i].end      = e;
                workers[i].has_work = 1;
                pthread_cond_signal(&workers[i].work_cv);
                pthread_mutex_unlock(&workers[i].mu);
            }
            for (int i = 0; i < N; i++) {
                pthread_mutex_lock(&workers[i].mu);
                while (workers[i].has_work) {
                    pthread_cond_wait(&workers[i].done_cv, &workers[i].mu);
                }
                pthread_mutex_unlock(&workers[i].mu);
            }
        }
        for (int i = 0; i < N; i++) {
            stats->adapter_trimmed_r1 += workers[i].local_trim_r1;
            stats->adapter_trimmed_r2 += workers[i].local_trim_r2;
        }

        int write_err = 0;
        for (int i = 0; i < n; i++) {
            if (account_and_emit(i, b->r1, b->r2, b->out, aw, fw, rw, stats, p) != 0) {
                fprintf(stderr, "pbj: write error\n");
                write_err = 1;
                break;
            }
        }

        pthread_mutex_lock(&ring.mu);
        ring.tail = (ring.tail + 1) % NSLOTS;
        ring.count--;
        if (write_err) {
            ring.fatal = 1;
            pthread_cond_broadcast(&ring.not_full_cv);
        } else {
            pthread_cond_signal(&ring.not_full_cv);
        }
        pthread_mutex_unlock(&ring.mu);
        if (write_err) { rc = 2; break; }
    }

done:
    if (reader_started) {
        pthread_mutex_lock(&ring.mu);
        ring.fatal = 1;
        pthread_cond_broadcast(&ring.not_full_cv);
        pthread_mutex_unlock(&ring.mu);
        pthread_join(reader_tid, NULL);
    }
    if (workers) {
        for (int i = 0; i < workers_started; i++) {
            pthread_mutex_lock(&workers[i].mu);
            workers[i].exit_flag = 1;
            pthread_cond_signal(&workers[i].work_cv);
            pthread_mutex_unlock(&workers[i].mu);
        }
        for (int i = 0; i < workers_started; i++) {
            pthread_join(workers[i].tid, NULL);
            pthread_cond_destroy(&workers[i].done_cv);
            pthread_cond_destroy(&workers[i].work_cv);
            pthread_mutex_destroy(&workers[i].mu);
        }
        for (int i = 0; i < N; i++) pbj_workspace_free(&workers[i].ws);
        free(workers);
    }
    if (ring_inited) {
        pthread_cond_destroy(&ring.not_empty_cv);
        pthread_cond_destroy(&ring.not_full_cv);
        pthread_mutex_destroy(&ring.mu);
    }
    for (int s = 0; s < NSLOTS; s++) {
        if (ring.slots[s].r1) {
            for (int i = 0; i < B; i++) pbj_read_free(&ring.slots[s].r1[i]);
            free(ring.slots[s].r1);
        }
        if (ring.slots[s].r2) {
            for (int i = 0; i < B; i++) pbj_read_free(&ring.slots[s].r2[i]);
            free(ring.slots[s].r2);
        }
        if (ring.slots[s].out) {
            for (int i = 0; i < B; i++) pbj_read_free(&ring.slots[s].out[i].merged);
            free(ring.slots[s].out);
        }
    }
    return rc;
}
