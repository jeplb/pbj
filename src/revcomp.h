#ifndef PBJ_REVCOMP_H
#define PBJ_REVCOMP_H

#include "pbj.h"

void pbj_revcomp_init(void);

/* writes the reverse complement of src->seq and reverse of src->qual into dst.
   dst is resized to fit. src and dst may not alias. */
void pbj_revcomp(const pbj_read_t *src, pbj_read_t *dst);

#endif
