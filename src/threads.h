#ifndef PBJ_THREADS_H
#define PBJ_THREADS_H

#include "pbj.h"
#include "fastq.h"

/* runs the merge pipeline single- or multi-threaded depending on p->threads.
   returns 0 on success, non-zero on error. final stats are written to *stats. */
int pbj_run_pipeline(const pbj_params_t *p,
                     pbj_fastq_reader_t *r1_reader,
                     pbj_fastq_reader_t *r2_reader,
                     pbj_fastq_writer_t *assembled_w,
                     pbj_fastq_writer_t *unassembled_fwd_w,
                     pbj_fastq_writer_t *unassembled_rev_w,
                     pbj_stats_t *stats);

#endif
