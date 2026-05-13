#ifndef PBJ_H
#define PBJ_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <zlib.h>

#define PBJ_VERSION "0.1.0"

/* result codes for pbj_try_merge */
typedef enum {
    PBJ_MERGED_OK         = 0,
    PBJ_NO_CANDIDATE      = 1,  /* no kmer match */
    PBJ_FAIL_PVALUE       = 2,
    PBJ_FAIL_LENGTH       = 3,
    PBJ_FAIL_AMBIGUOUS    = 4,
    PBJ_FAIL_MIN_QUAL     = 5,
    PBJ_FAIL_MAX_EE       = 6,  /* expected errors above --max-ee */
    PBJ_FAIL_MONO         = 7   /* low-complexity gate, single base too dominant */
} pbj_status_t;

/* cli parameters and tuning knobs */
typedef struct {
    /* i/o */
    const char *r1_path;
    const char *r2_path;
    const char *out_prefix;
    int         gzip_output;
    int         gzip_level;
    int         interleaved_in;     /* read alternating pairs from r1_path */

    /* algorithm */
    int    min_overlap;
    int    max_overlap;
    int    min_assembly_len;
    int    max_assembly_len;
    int    kmer_size;
    int    match_reward;
    int    mismatch_penalty;
    double p_value;
    double p_match;
    int    bonferroni;
    int    quality_offset;
    int    min_mean_qual_merged;
    int    strict_ambiguity;
    int    strict_q_threshold;
    int    allow_dovetail;

    /* adapter trimming, applied before merging */
    const char *adapter_forward;
    const char *adapter_reverse;
    int    adapter_min_match;        /* min bp of adapter overlap */
    double adapter_match_frac;       /* min match fraction over the overlap */

    /* post-merge filters */
    double max_ee;                   /* 0 = disabled; max expected errors */
    double max_mono_frac;            /* 0 = disabled; max single-base fraction */
    int    qtrim;                    /* 0 = disabled; trim 3' bases below phred n */

    /* runtime */
    int    threads;
    int    batch_size;
    int    verbose;
} pbj_params_t;

typedef struct {
    char     *name;
    uint32_t  name_len;
    char     *seq;
    char     *qual;
    uint32_t  len;
    uint32_t  cap;
} pbj_read_t;

typedef struct {
    pbj_status_t status;
    pbj_read_t   merged;
    int          offset;
    int          overlap_len;
    int          matches;
    int          mismatches;
    int          n_candidates;
    double       score;
    double       p_value;
    double       expected_errors;    /* valid when status == pbj_merged_ok */
} pbj_merge_result_t;

typedef struct {
    pbj_read_t r2_rc;
    void      *kmer_index;
    int       *candidates;
    int        candidates_cap;
    uint8_t   *cand_bitset;
    int        cand_bitset_bytes;
} pbj_workspace_t;

typedef struct {
    uint64_t total_pairs;
    uint64_t merged;
    uint64_t no_candidate;
    uint64_t fail_pvalue;
    uint64_t fail_length;
    uint64_t fail_ambiguous;
    uint64_t fail_min_qual;
    uint64_t fail_max_ee;
    uint64_t fail_mono;
    uint64_t adapter_trimmed_r1;
    uint64_t adapter_trimmed_r2;
} pbj_stats_t;

void pbj_params_init_defaults(pbj_params_t *p);
void pbj_workspace_init(pbj_workspace_t *ws, const pbj_params_t *p);
void pbj_workspace_free(pbj_workspace_t *ws);
void pbj_read_init(pbj_read_t *r);
void pbj_read_free(pbj_read_t *r);
void pbj_read_ensure_cap(pbj_read_t *r, uint32_t cap);

int pbj_pair_names_match(const pbj_read_t *r1, const pbj_read_t *r2);

/* returns 1 if the read's casava 1.8+ comment begins with "1:" (i.e., it
   is the forward member of a pair), 2 if it begins with "2:", or 0 if no
   recognizable directional flag is present. */
int pbj_casava_direction(const pbj_read_t *r);

#endif
