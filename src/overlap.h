#ifndef PBJ_OVERLAP_H
#define PBJ_OVERLAP_H

#include "pbj.h"

void pbj_overlap_init(pbj_workspace_t *ws, const pbj_params_t *p);
void pbj_overlap_free(pbj_workspace_t *ws);

/* generates candidate signed offsets s where r2_rc[max(0,-s) ...] aligns
   with r1[max(0,s) ...]. results are deduplicated and stored in
   ws->candidates. returns the count. */
int pbj_overlap_candidates(pbj_workspace_t *ws,
                            const pbj_read_t *r1,
                            const pbj_read_t *r2_rc,
                            const pbj_params_t *p);

/* given offset s, computes the overlap geometry:
     r1_start, r2_start, overlap_len
   returns 1 if the overlap is well-formed and satisfies min/max_overlap. */
int pbj_overlap_geometry(int s,
                          int len_r1, int len_r2,
                          int min_overlap, int max_overlap,
                          int *r1_start, int *r2_start, int *olen);

#endif
