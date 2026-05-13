#ifndef PBJ_MERGE_H
#define PBJ_MERGE_H

#include "pbj.h"

/* tries to merge a pair. fills result and, if successful, the merged read. */
void pbj_try_merge(const pbj_read_t *r1,
                   const pbj_read_t *r2,
                   pbj_workspace_t *ws,
                   const pbj_params_t *p,
                   pbj_merge_result_t *result);

#endif
