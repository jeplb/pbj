#ifndef PBJ_FASTQ_H
#define PBJ_FASTQ_H

#include "pbj.h"

typedef struct pbj_fastq_reader pbj_fastq_reader_t;
typedef struct pbj_fastq_writer pbj_fastq_writer_t;

pbj_fastq_reader_t *pbj_fastq_open_reader(const char *path);
void                pbj_fastq_close_reader(pbj_fastq_reader_t *r);

/* returns 1 on success, 0 on eof, -1 on error */
int pbj_fastq_read(pbj_fastq_reader_t *r, pbj_read_t *out);

/* gzip_level is 1..9 when use_gzip is on, ignored otherwise. */
pbj_fastq_writer_t *pbj_fastq_open_writer(const char *path, int use_gzip, int gzip_level);

/* closes and flushes the writer. returns 0 on success, -1 if the final
   flush failed (disk full, broken pipe, gzip dictionary error). */
int  pbj_fastq_close_writer(pbj_fastq_writer_t *w);

/* writes a single record. returns 0 on success, -1 on error. output
   quality bytes are translated from input_offset to Phred+33. */
int pbj_fastq_write(pbj_fastq_writer_t *w, const pbj_read_t *rec, int input_offset);

#endif
