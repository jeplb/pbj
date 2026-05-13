#ifndef PBJ_ADAPTER_H
#define PBJ_ADAPTER_H

#include "pbj.h"

/* find the leftmost position in seq where a suffix-aligned adapter match
   reaches the required quality. returns len if no acceptable match exists.
   the returned position is the new (post-trim) read length. */
int pbj_adapter_find_trim(const char *seq, int len,
                          const char *adapter, int adapter_len,
                          int min_match_bp, double min_match_frac);

/* in-place trim. returns 1 if a trim occurred, 0 otherwise. */
int pbj_adapter_trim(pbj_read_t *r,
                     const char *adapter,
                     int min_match_bp, double min_match_frac);

#endif
