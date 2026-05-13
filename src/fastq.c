#include "fastq.h"
#include "kseq.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

struct pbj_fastq_reader {
    gzFile  fp;
    kseq_t *ks;
};

struct pbj_fastq_writer {
    int    use_gzip;
    gzFile gz;
    FILE  *fp;
};

void pbj_read_init(pbj_read_t *r) {
    r->name     = NULL;
    r->name_len = 0;
    r->seq      = NULL;
    r->qual     = NULL;
    r->len      = 0;
    r->cap      = 0;
}

void pbj_read_free(pbj_read_t *r) {
    free(r->name);
    free(r->seq);
    free(r->qual);
    r->name     = NULL;
    r->name_len = 0;
    r->seq      = NULL;
    r->qual     = NULL;
    r->len      = 0;
    r->cap      = 0;
}

void pbj_read_ensure_cap(pbj_read_t *r, uint32_t cap) {
    if (r->cap >= cap) return;
    uint32_t newcap = r->cap ? r->cap : 64;
    while (newcap < cap) newcap *= 2;
    /* commit each pointer immediately after its realloc succeeds. on second
       failure the first new buffer is reachable via r->seq (not dangling)
       and r->qual still points at the unfreed old buffer. avoids a uaf
       window if abort is ever swapped for a recovery path. */
    char *new_seq = (char*)realloc(r->seq, newcap);
    if (!new_seq) {
        fprintf(stderr, "pbj: out of memory in pbj_read_ensure_cap (%u bytes)\n", newcap);
        abort();
    }
    r->seq = new_seq;
    char *new_qual = (char*)realloc(r->qual, newcap);
    if (!new_qual) {
        fprintf(stderr, "pbj: out of memory in pbj_read_ensure_cap (%u bytes)\n", newcap);
        abort();
    }
    r->qual = new_qual;
    r->cap  = newcap;
}

pbj_fastq_reader_t *pbj_fastq_open_reader(const char *path) {
    gzFile fp;
    if (path && (strcmp(path, "-") == 0)) {
        fp = gzdopen(fileno(stdin), "r");
    } else {
        fp = gzopen(path, "r");
    }
    if (!fp) return NULL;
    pbj_fastq_reader_t *r = (pbj_fastq_reader_t*)calloc(1, sizeof(*r));
    if (!r) { gzclose(fp); return NULL; }
    r->fp = fp;
    r->ks = kseq_init(fp);
    return r;
}

void pbj_fastq_close_reader(pbj_fastq_reader_t *r) {
    if (!r) return;
    if (r->ks) kseq_destroy(r->ks);
    if (r->fp) gzclose(r->fp);
    free(r);
}

/* hard cap on a single record length to keep all int-arithmetic in
   overlap.c / merge.c well within int_max and to avoid pathological
   memory use. realistic paired-end reads never approach this. */
#define PBJ_MAX_READ_LEN (1 << 28)  /* 256 mib */

int pbj_fastq_read(pbj_fastq_reader_t *r, pbj_read_t *out) {
    int n = kseq_read(r->ks);
    if (n < 0) {
        if (n == -1) return 0;
        return -1;
    }
    if (n > PBJ_MAX_READ_LEN) {
        fprintf(stderr, "pbj: read exceeds maximum supported length (%d > %d)\n",
                n, PBJ_MAX_READ_LEN);
        return -1;
    }
    size_t name_len    = r->ks->name.l;
    size_t comment_len = r->ks->comment.l;
    size_t total_name  = name_len + (comment_len ? comment_len + 1 : 0);
    char *new_name = (char*)realloc(out->name, total_name + 1);
    if (!new_name) return -1;
    out->name = new_name;
    memcpy(out->name, r->ks->name.s, name_len);
    if (comment_len) {
        out->name[name_len] = ' ';
        memcpy(out->name + name_len + 1, r->ks->comment.s, comment_len);
    }
    out->name[total_name] = '\0';
    out->name_len = (uint32_t)name_len;

    pbj_read_ensure_cap(out, (uint32_t)n + 1);
    /* case-normalize seq to uppercase; downstream comparisons assume this */
    for (int i = 0; i < n; i++) {
        char c = r->ks->seq.s[i];
        if (c >= 'a' && c <= 'z') c = (char)(c - 32);
        out->seq[i] = c;
    }
    out->seq[n] = '\0';
    if (r->ks->qual.l == (size_t)n) {
        memcpy(out->qual, r->ks->qual.s, n);
    } else {
        /* fasta input, synthesize quality */
        memset(out->qual, '!' + 30, n);
    }
    out->qual[n] = '\0';
    out->len = (uint32_t)n;
    return 1;
}

/* core-id comparison: name up to first space (or name_len), with optional
   /1 or /2 mate-suffix stripped. tracks the stripped suffix character so
   that callers can reject pairs that stripped the same explicit suffix
   (i.e., r1 and r2 both ending in /1, or both in /2). */
typedef struct { uint32_t len; char stripped; } core_info_t;

static core_info_t core_info(const pbj_read_t *r) {
    core_info_t out = { r->name_len, 0 };
    if (r->name_len >= 2 && r->name[r->name_len - 2] == '/' &&
        (r->name[r->name_len - 1] == '1' || r->name[r->name_len - 1] == '2')) {
        out.stripped = r->name[r->name_len - 1];
        out.len     = r->name_len - 2;
    }
    return out;
}

int pbj_pair_names_match(const pbj_read_t *r1, const pbj_read_t *r2) {
    core_info_t a = core_info(r1);
    core_info_t b = core_info(r2);
    if (a.len != b.len) return 0;
    if (a.len > 0 && memcmp(r1->name, r2->name, a.len) != 0) return 0;
    /* both stripped the same explicit suffix means these are not a pair
       (e.g., both records ended with /1). */
    if (a.stripped != 0 && a.stripped == b.stripped) return 0;
    return 1;
}

int pbj_casava_direction(const pbj_read_t *r) {
    if (!r->name) return 0;
    /* the comment field begins one byte after name_len when present;
       fastq.c stores name and comment joined by a space. */
    uint32_t n = (uint32_t)strlen(r->name);
    if (n <= r->name_len + 2)   return 0;
    if (r->name[r->name_len] != ' ') return 0;
    char d = r->name[r->name_len + 1];
    if ((d != '1' && d != '2')) return 0;
    if (r->name[r->name_len + 2] != ':') return 0;
    return d - '0';
}

pbj_fastq_writer_t *pbj_fastq_open_writer(const char *path, int use_gzip, int gzip_level) {
    pbj_fastq_writer_t *w = (pbj_fastq_writer_t*)calloc(1, sizeof(*w));
    if (!w) return NULL;
    w->use_gzip = use_gzip;
    if (use_gzip) {
        if (gzip_level < 1) gzip_level = 1;
        if (gzip_level > 9) gzip_level = 9;
        char mode[8];
        snprintf(mode, sizeof(mode), "wb%d", gzip_level);
        w->gz = gzopen(path, mode);
        if (!w->gz) { free(w); return NULL; }
    } else {
        w->fp = fopen(path, "wb");
        if (!w->fp) { free(w); return NULL; }
        setvbuf(w->fp, NULL, _IOFBF, 1 << 20);
    }
    return w;
}

int pbj_fastq_close_writer(pbj_fastq_writer_t *w) {
    if (!w) return 0;
    int rc = 0;
    if (w->gz) {
        if (gzclose(w->gz) != Z_OK) rc = -1;
    }
    if (w->fp) {
        if (fflush(w->fp) != 0)    rc = -1;
        if (fclose(w->fp) != 0)    rc = -1;
    }
    free(w);
    return rc;
}

/* writes len qual bytes translated from input_offset to phred+33. */
static int write_qual_translated(pbj_fastq_writer_t *w,
                                  const char *qual, uint32_t len, int input_offset) {
    char buf[4096];
    int delta = 33 - input_offset;
    if (delta == 0) {
        if (w->use_gzip) return gzwrite(w->gz, qual, len) == (int)len ? 0 : -1;
        return fwrite(qual, 1, len, w->fp) == len ? 0 : -1;
    }
    uint32_t i = 0;
    while (i < len) {
        uint32_t chunk = (len - i) < sizeof(buf) ? (len - i) : (uint32_t)sizeof(buf);
        for (uint32_t j = 0; j < chunk; j++) {
            int q = (int)(unsigned char)qual[i + j] + delta;
            if (q < 33)  q = 33;
            if (q > 126) q = 126;
            buf[j] = (char)q;
        }
        if (w->use_gzip) {
            if (gzwrite(w->gz, buf, chunk) != (int)chunk) return -1;
        } else {
            if (fwrite(buf, 1, chunk, w->fp) != chunk) return -1;
        }
        i += chunk;
    }
    return 0;
}

int pbj_fastq_write(pbj_fastq_writer_t *w, const pbj_read_t *rec, int input_offset) {
    if (w->use_gzip) {
        if (gzputc(w->gz, '@') < 0) return -1;
        if (gzputs(w->gz, rec->name ? rec->name : "read") < 0) return -1;
        if (gzputc(w->gz, '\n') < 0) return -1;
        if (gzwrite(w->gz, rec->seq, rec->len) != (int)rec->len) return -1;
        if (gzwrite(w->gz, "\n+\n", 3) != 3) return -1;
        if (write_qual_translated(w, rec->qual, rec->len, input_offset) != 0) return -1;
        if (gzputc(w->gz, '\n') < 0) return -1;
    } else {
        if (fputc('@', w->fp) == EOF) return -1;
        if (fputs(rec->name ? rec->name : "read", w->fp) == EOF) return -1;
        if (fputc('\n', w->fp) == EOF) return -1;
        if (fwrite(rec->seq, 1, rec->len, w->fp) != rec->len) return -1;
        if (fwrite("\n+\n", 1, 3, w->fp) != 3) return -1;
        if (write_qual_translated(w, rec->qual, rec->len, input_offset) != 0) return -1;
        if (fputc('\n', w->fp) == EOF) return -1;
    }
    return 0;
}
