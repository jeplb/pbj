#include "pbj.h"
#include "fastq.h"
#include "revcomp.h"
#include "overlap.h"
#include "threads.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <errno.h>

static void usage(FILE *out) {
    fprintf(out,
        "pbj %s, paired-end read merger\n"
        "usage: pbj -f R1.fq[.gz] [-r R2.fq[.gz] | --interleaved-in] -o PREFIX [options]\n\n"
        "required:\n"
        "  -f, --forward FILE        forward reads (R1), or interleaved input; '-' = stdin\n"
        "  -r, --reverse FILE        reverse reads (R2), omit when --interleaved-in; '-' = stdin\n"
        "  -o, --output PREFIX       output prefix\n\n"
        "io:\n"
        "  -z, --gzip                gzip-compress output\n"
        "      --gzip-level N        gzip level 1..9 [6]\n"
        "      --interleaved-in      read R1 and R2 alternating from -f\n\n"
        "algorithm:\n"
        "  -m, --min-overlap N       minimum overlap length [20]\n"
        "  -M, --max-overlap N       maximum overlap length [0 = unbounded]\n"
        "  -l, --min-length N        minimum merged length [50]\n"
        "  -L, --max-length N        maximum merged length [0 = unbounded]\n"
        "  -k, --kmer-size N         kmer size for seeding [12]\n"
        "  -p, --p-value F           p-value cutoff for accepting a merge [0.01]\n"
        "      --p-match F           null match probability [0.25]\n"
        "      --no-bonferroni       skip Bonferroni correction over candidates\n"
        "  -q, --quality-offset N    input Phred offset, 33 or 64 [33]; output is +33\n"
        "      --match REWARD        match reward for scoring [1]\n"
        "      --mismatch PENALTY    mismatch penalty for scoring [1]\n"
        "  -s, --strict              reject overlaps with high-Q mismatches\n"
        "      --strict-q N          Phred above which a mismatch is high-Q [20]\n"
        "  -d, --dovetail            allow dovetail / adapter-readthrough alignments\n\n"
        "filters (applied to merged read):\n"
        "  -Q, --min-mean-qual N     minimum mean Phred of merged read [0 = off]\n"
        "      --max-ee F            max expected errors per merged read [0 = off]\n"
        "      --max-mono-frac F     max single-base fraction in [0,1] [0 = off]\n"
        "      --qtrim N             trim 3' bases below Phred N [0 = off]\n\n"
        "adapter trimming (applied to input reads before merge):\n"
        "  -a, --adapter-forward SEQ adapter expected at R1 3' end\n"
        "  -A, --adapter-reverse SEQ adapter expected at R2 3' end\n"
        "      --adapter-min-match N min bp of adapter overlap [8]\n"
        "      --adapter-match-frac F min match fraction [0.9]\n\n"
        "runtime:\n"
        "  -t, --threads N           worker threads [1]\n"
        "  -b, --batch-size N        pairs per work batch [1024]\n"
        "  -v, --verbose             per-pair status to stderr\n"
        "  -h, --help                this help\n"
        "      --version             version\n",
        PBJ_VERSION);
}

void pbj_params_init_defaults(pbj_params_t *p) {
    memset(p, 0, sizeof(*p));
    p->min_overlap          = 20;
    p->max_overlap          = 0;
    p->min_assembly_len     = 50;
    p->max_assembly_len     = 0;
    p->kmer_size            = 12;
    p->match_reward         = 1;
    p->mismatch_penalty     = 1;
    p->p_value              = 0.01;
    p->p_match              = 0.25;
    p->bonferroni           = 1;
    p->quality_offset       = 33;
    p->min_mean_qual_merged = 0;
    p->strict_ambiguity     = 0;
    p->strict_q_threshold   = 20;
    p->allow_dovetail       = 0;
    p->adapter_forward      = NULL;
    p->adapter_reverse      = NULL;
    p->adapter_min_match    = 8;
    p->adapter_match_frac   = 0.9;
    p->max_ee               = 0.0;
    p->max_mono_frac        = 0.0;
    p->qtrim                = 0;
    p->threads              = 1;
    p->batch_size           = 1024;
    p->verbose              = 0;
    p->gzip_output          = 0;
    p->gzip_level           = 6;
    p->interleaved_in       = 0;
}

void pbj_workspace_init(pbj_workspace_t *ws, const pbj_params_t *p) {
    memset(ws, 0, sizeof(*ws));
    pbj_read_init(&ws->r2_rc);
    pbj_overlap_init(ws, p);
}

void pbj_workspace_free(pbj_workspace_t *ws) {
    pbj_read_free(&ws->r2_rc);
    pbj_overlap_free(ws);
}

int main(int argc, char **argv) {
    pbj_params_t params;
    pbj_params_init_defaults(&params);

    static struct option long_opts[] = {
        {"forward",            required_argument, 0, 'f'},
        {"reverse",            required_argument, 0, 'r'},
        {"output",             required_argument, 0, 'o'},
        {"gzip",               no_argument,       0, 'z'},
        {"gzip-level",         required_argument, 0, 1007},
        {"interleaved-in",     no_argument,       0, 1008},
        {"min-overlap",        required_argument, 0, 'm'},
        {"max-overlap",        required_argument, 0, 'M'},
        {"min-length",         required_argument, 0, 'l'},
        {"max-length",         required_argument, 0, 'L'},
        {"kmer-size",          required_argument, 0, 'k'},
        {"p-value",            required_argument, 0, 'p'},
        {"p-match",            required_argument, 0, 1004},
        {"no-bonferroni",      no_argument,       0, 1005},
        {"quality-offset",     required_argument, 0, 'q'},
        {"min-mean-qual",      required_argument, 0, 'Q'},
        {"match",              required_argument, 0, 1001},
        {"mismatch",           required_argument, 0, 1002},
        {"strict",             no_argument,       0, 's'},
        {"strict-q",           required_argument, 0, 1006},
        {"dovetail",           no_argument,       0, 'd'},
        {"adapter-forward",    required_argument, 0, 'a'},
        {"adapter-reverse",    required_argument, 0, 'A'},
        {"adapter-min-match",  required_argument, 0, 1009},
        {"adapter-match-frac", required_argument, 0, 1010},
        {"max-ee",             required_argument, 0, 1011},
        {"max-mono-frac",      required_argument, 0, 1012},
        {"qtrim",              required_argument, 0, 1013},
        {"threads",            required_argument, 0, 't'},
        {"batch-size",         required_argument, 0, 'b'},
        {"verbose",            no_argument,       0, 'v'},
        {"help",               no_argument,       0, 'h'},
        {"version",            no_argument,       0, 1003},
        {0,0,0,0}
    };

    int c;
    while ((c = getopt_long(argc, argv,
                            "f:r:o:zm:M:l:L:k:p:q:Q:sda:A:t:b:vh",
                            long_opts, NULL)) != -1) {
        switch (c) {
            case 'f': params.r1_path             = optarg;        break;
            case 'r': params.r2_path             = optarg;        break;
            case 'o': params.out_prefix          = optarg;        break;
            case 'z': params.gzip_output         = 1;             break;
            case 'm': params.min_overlap         = atoi(optarg);  break;
            case 'M': params.max_overlap         = atoi(optarg);  break;
            case 'l': params.min_assembly_len    = atoi(optarg);  break;
            case 'L': params.max_assembly_len    = atoi(optarg);  break;
            case 'k': params.kmer_size           = atoi(optarg);  break;
            case 'p': params.p_value             = atof(optarg);  break;
            case 'q': params.quality_offset      = atoi(optarg);  break;
            case 'Q': params.min_mean_qual_merged= atoi(optarg);  break;
            case 'a': params.adapter_forward     = optarg;        break;
            case 'A': params.adapter_reverse     = optarg;        break;
            case 1001: params.match_reward       = atoi(optarg);  break;
            case 1002: params.mismatch_penalty   = atoi(optarg);  break;
            case 1003: printf("pbj %s\n", PBJ_VERSION); return 0;
            case 1004: params.p_match            = atof(optarg);  break;
            case 1005: params.bonferroni         = 0;             break;
            case 1006: params.strict_q_threshold = atoi(optarg);  break;
            case 1007: params.gzip_level         = atoi(optarg);  break;
            case 1008: params.interleaved_in     = 1;             break;
            case 1009: params.adapter_min_match  = atoi(optarg);  break;
            case 1010: params.adapter_match_frac = atof(optarg);  break;
            case 1011: params.max_ee             = atof(optarg);  break;
            case 1012: params.max_mono_frac      = atof(optarg);  break;
            case 1013: params.qtrim              = atoi(optarg);  break;
            case 's': params.strict_ambiguity    = 1;             break;
            case 'd': params.allow_dovetail      = 1;             break;
            case 't': params.threads             = atoi(optarg);  break;
            case 'b': params.batch_size          = atoi(optarg);  break;
            case 'v': params.verbose             = 1;             break;
            case 'h': usage(stdout); return 0;
            default: usage(stderr); return 1;
        }
    }

    if (!params.r1_path || !params.out_prefix) {
        fprintf(stderr, "error: -f and -o are required\n\n");
        usage(stderr);
        return 1;
    }
    if (params.interleaved_in) {
        if (params.r2_path) {
            fprintf(stderr, "error: --interleaved-in is incompatible with -r\n");
            return 1;
        }
    } else if (!params.r2_path) {
        fprintf(stderr, "error: -r is required unless --interleaved-in is set\n");
        return 1;
    }
    if (params.threads < 1)        params.threads    = 1;
    if (params.batch_size < 1)     params.batch_size = 1;
    if (params.kmer_size < 4 || params.kmer_size > 31) {
        fprintf(stderr, "error: kmer-size must be in [4, 31]\n"); return 1;
    }
    if (params.p_value < 0.0 || params.p_value > 1.0) {
        fprintf(stderr, "error: p-value must be in [0, 1]\n"); return 1;
    }
    if (params.p_match <= 0.0 || params.p_match >= 1.0) {
        fprintf(stderr, "error: p-match must be in (0, 1)\n"); return 1;
    }
    if (params.gzip_level < 1 || params.gzip_level > 9) {
        fprintf(stderr, "error: gzip-level must be in [1, 9]\n"); return 1;
    }
    if (params.max_ee < 0.0) {
        fprintf(stderr, "error: max-ee must be >= 0\n"); return 1;
    }
    if (params.max_mono_frac < 0.0 || params.max_mono_frac > 1.0) {
        fprintf(stderr, "error: max-mono-frac must be in [0, 1]\n"); return 1;
    }
    if (params.qtrim < 0) {
        fprintf(stderr, "error: qtrim must be >= 0\n"); return 1;
    }
    if (params.quality_offset != 33 && params.quality_offset != 64) {
        fprintf(stderr, "warning: unusual quality offset %d (expected 33 or 64)\n",
                params.quality_offset);
    }
    if (params.adapter_forward && !params.adapter_reverse) {
        params.adapter_reverse = params.adapter_forward;
    } else if (params.adapter_reverse && !params.adapter_forward) {
        params.adapter_forward = params.adapter_reverse;
    }

    pbj_revcomp_init();

    pbj_fastq_reader_t *r1 = pbj_fastq_open_reader(params.r1_path);
    if (!r1) { fprintf(stderr, "error: cannot open %s: %s\n", params.r1_path, strerror(errno)); return 2; }
    pbj_fastq_reader_t *r2;
    if (params.interleaved_in) {
        /* alternating reads come from r1 */
        r2 = r1;
    } else {
        r2 = pbj_fastq_open_reader(params.r2_path);
        if (!r2) {
            fprintf(stderr, "error: cannot open %s: %s\n", params.r2_path, strerror(errno));
            pbj_fastq_close_reader(r1);
            return 2;
        }
    }

    char path[2048];
    const char *ext = params.gzip_output ? ".fastq.gz" : ".fastq";

    /* check each snprintf for truncation so a very long --output prefix
       can't silently produce colliding output paths. */
    int wlen = snprintf(path, sizeof(path), "%s.assembled%s", params.out_prefix, ext);
    if (wlen < 0 || (size_t)wlen >= sizeof(path)) {
        fprintf(stderr, "error: --output prefix is too long\n");
        pbj_fastq_close_reader(r1);
        if (r2 != r1) pbj_fastq_close_reader(r2);
        return 1;
    }
    pbj_fastq_writer_t *aw = pbj_fastq_open_writer(path, params.gzip_output, params.gzip_level);
    if (!aw) {
        fprintf(stderr, "error: cannot open %s\n", path);
        pbj_fastq_close_reader(r1);
        if (r2 != r1) pbj_fastq_close_reader(r2);
        return 2;
    }

    wlen = snprintf(path, sizeof(path), "%s.unassembled.forward%s", params.out_prefix, ext);
    if (wlen < 0 || (size_t)wlen >= sizeof(path)) {
        fprintf(stderr, "error: --output prefix is too long\n");
        pbj_fastq_close_writer(aw); pbj_fastq_close_reader(r1);
        if (r2 != r1) pbj_fastq_close_reader(r2);
        return 1;
    }
    pbj_fastq_writer_t *fw = pbj_fastq_open_writer(path, params.gzip_output, params.gzip_level);
    if (!fw) {
        fprintf(stderr, "error: cannot open %s\n", path);
        pbj_fastq_close_writer(aw); pbj_fastq_close_reader(r1);
        if (r2 != r1) pbj_fastq_close_reader(r2);
        return 2;
    }

    wlen = snprintf(path, sizeof(path), "%s.unassembled.reverse%s", params.out_prefix, ext);
    if (wlen < 0 || (size_t)wlen >= sizeof(path)) {
        fprintf(stderr, "error: --output prefix is too long\n");
        pbj_fastq_close_writer(aw); pbj_fastq_close_writer(fw); pbj_fastq_close_reader(r1);
        if (r2 != r1) pbj_fastq_close_reader(r2);
        return 1;
    }
    pbj_fastq_writer_t *rw = pbj_fastq_open_writer(path, params.gzip_output, params.gzip_level);
    if (!rw) {
        fprintf(stderr, "error: cannot open %s\n", path);
        pbj_fastq_close_writer(aw); pbj_fastq_close_writer(fw); pbj_fastq_close_reader(r1);
        if (r2 != r1) pbj_fastq_close_reader(r2);
        return 2;
    }

    pbj_stats_t stats;
    memset(&stats, 0, sizeof(stats));

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int rc = pbj_run_pipeline(&params, r1, r2, aw, fw, rw, &stats);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double seconds = (double)(t1.tv_sec - t0.tv_sec) + 1e-9 * (double)(t1.tv_nsec - t0.tv_nsec);

    pbj_fastq_close_reader(r1);
    if (r2 != r1) pbj_fastq_close_reader(r2);
    int wa = pbj_fastq_close_writer(aw);
    int wf = pbj_fastq_close_writer(fw);
    int wr = pbj_fastq_close_writer(rw);
    if (wa || wf || wr) {
        fprintf(stderr, "pbj: write error during flush/close (disk full or broken pipe)\n");
        if (rc == 0) rc = 2;
    }

    fprintf(stderr, "pbj summary:\n");
    fprintf(stderr, "  total pairs        %12llu\n",
            (unsigned long long)stats.total_pairs);
    if (stats.total_pairs > 0) {
        double tp = (double)stats.total_pairs;
        fprintf(stderr, "  merged             %12llu (%6.2f%%)\n",
                (unsigned long long)stats.merged,        100.0 * (double)stats.merged        / tp);
        fprintf(stderr, "  no candidate       %12llu (%6.2f%%)\n",
                (unsigned long long)stats.no_candidate,  100.0 * (double)stats.no_candidate  / tp);
        fprintf(stderr, "  fail p-value       %12llu (%6.2f%%)\n",
                (unsigned long long)stats.fail_pvalue,   100.0 * (double)stats.fail_pvalue   / tp);
        fprintf(stderr, "  fail length        %12llu (%6.2f%%)\n",
                (unsigned long long)stats.fail_length,   100.0 * (double)stats.fail_length   / tp);
        fprintf(stderr, "  fail ambiguity     %12llu (%6.2f%%)\n",
                (unsigned long long)stats.fail_ambiguous,100.0 * (double)stats.fail_ambiguous/ tp);
        fprintf(stderr, "  fail min quality   %12llu (%6.2f%%)\n",
                (unsigned long long)stats.fail_min_qual, 100.0 * (double)stats.fail_min_qual / tp);
        fprintf(stderr, "  fail max-ee        %12llu (%6.2f%%)\n",
                (unsigned long long)stats.fail_max_ee,   100.0 * (double)stats.fail_max_ee   / tp);
        fprintf(stderr, "  fail low-complex   %12llu (%6.2f%%)\n",
                (unsigned long long)stats.fail_mono,     100.0 * (double)stats.fail_mono     / tp);
        if (params.adapter_forward || params.adapter_reverse) {
            fprintf(stderr, "  adapter trimmed R1 %12llu\n",
                    (unsigned long long)stats.adapter_trimmed_r1);
            fprintf(stderr, "  adapter trimmed R2 %12llu\n",
                    (unsigned long long)stats.adapter_trimmed_r2);
        }
        fprintf(stderr, "  throughput         %12.0f pairs/s\n",
                seconds > 0 ? tp / seconds : 0.0);
    }
    fprintf(stderr, "  elapsed            %12.3f s\n", seconds);

    return rc;
}
