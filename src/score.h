#ifndef PBJ_SCORE_H
#define PBJ_SCORE_H

#include "pbj.h"

/* quality-weighted score plus raw match/mismatch counts for the overlap
   region. caller has already shifted offsets via pbj_overlap_geometry. */
typedef struct {
    int    matches;
    int    mismatches;
    int    n_count;     /* positions where either base is n */
    double score;
} pbj_score_t;

void pbj_score_overlap(const pbj_read_t *r1, int r1_start,
                        const pbj_read_t *r2_rc, int r2_start,
                        int overlap_len,
                        const pbj_params_t *p,
                        pbj_score_t *out);

#endif
