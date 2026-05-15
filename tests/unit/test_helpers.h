/* shared assertion macros and helper functions for pbj unit-level tests.
   each including .c file gets its own copy of g_failed / g_tests because
   they're declared static here; this keeps the file as a single-binary
   helper without needing a separate .c. */

#ifndef PBJ_TEST_HELPERS_H
#define PBJ_TEST_HELPERS_H

#include "../../src/pbj.h"
#include "../../src/revcomp.h"
#include "../../src/stats.h"
#include "../../src/score.h"
#include "../../src/overlap.h"
#include "../../src/adapter.h"
#include "../../src/merge.h"
#include "../../src/fastq.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int g_failed = 0;
static int g_tests  = 0;
static const char *g_current = "";

#define CHECK(cond) do {                                                  \
    g_tests++;                                                             \
    if (!(cond)) {                                                         \
        g_failed++;                                                        \
        fprintf(stderr, "  FAIL %s: %s:%d: %s\n",                          \
                g_current, __FILE__, __LINE__, #cond);                     \
    }                                                                      \
} while (0)

#define CHECK_EQ_INT(a, b) do {                                            \
    g_tests++;                                                             \
    long _a = (long)(a), _b = (long)(b);                                   \
    if (_a != _b) {                                                        \
        g_failed++;                                                        \
        fprintf(stderr, "  FAIL %s: %s:%d: %s (%ld) != %s (%ld)\n",        \
                g_current, __FILE__, __LINE__, #a, _a, #b, _b);            \
    }                                                                      \
} while (0)

#define CHECK_EQ_STR(a, b) do {                                            \
    g_tests++;                                                             \
    const char *_a = (a), *_b = (b);                                       \
    if (strcmp(_a, _b) != 0) {                                             \
        g_failed++;                                                        \
        fprintf(stderr, "  FAIL %s: %s:%d: %s ('%s') != %s ('%s')\n",      \
                g_current, __FILE__, __LINE__, #a, _a, #b, _b);            \
    }                                                                      \
} while (0)

#define CHECK_CLOSE(a, b, eps) do {                                        \
    g_tests++;                                                             \
    double _a = (double)(a), _b = (double)(b);                             \
    if (fabs(_a - _b) > (eps)) {                                           \
        g_failed++;                                                        \
        fprintf(stderr, "  FAIL %s: %s:%d: |%s - %s| = %g > %g\n",         \
                g_current, __FILE__, __LINE__, #a, #b,                     \
                fabs(_a - _b), (double)(eps));                             \
    }                                                                      \
} while (0)

#define RUN(fn) do { g_current = #fn; fn(); } while (0)

static void params_defaults(pbj_params_t *p) {
    memset(p, 0, sizeof(*p));
    p->min_overlap          = 20;
    p->min_assembly_len     = 50;
    p->kmer_size            = 12;
    p->match_reward         = 1;
    p->mismatch_penalty     = 1;
    p->p_value              = 0.01;
    p->p_match              = 0.25;
    p->bonferroni           = 1;
    p->quality_offset       = 33;
    p->strict_q_threshold   = 20;
    p->adapter_min_match    = 8;
    p->adapter_match_frac   = 0.9;
    p->gzip_level           = 6;
    p->threads              = 1;
    p->batch_size           = 1024;
}

static void ws_init(pbj_workspace_t *ws, const pbj_params_t *p) {
    memset(ws, 0, sizeof(*ws));
    pbj_read_init(&ws->r2_rc);
    pbj_overlap_init(ws, p);
}

static void ws_free(pbj_workspace_t *ws) {
    pbj_read_free(&ws->r2_rc);
    pbj_overlap_free(ws);
}

static void load_read(pbj_read_t *r, const char *name,
                      const char *seq, const char *qual) {
    pbj_read_init(r);
    uint32_t L = (uint32_t)strlen(seq);
    pbj_read_ensure_cap(r, L + 1);
    memcpy(r->seq, seq, L); r->seq[L] = '\0';
    if (qual) {
        memcpy(r->qual, qual, L); r->qual[L] = '\0';
    } else {
        memset(r->qual, '!' + 30, L); r->qual[L] = '\0';
    }
    r->len = L;
    size_t nlen = strlen(name);
    r->name = (char*)malloc(nlen + 1);
    memcpy(r->name, name, nlen + 1);
    r->name_len = (uint32_t)nlen;
}

static char *revcomp_str(const char *s) {
    size_t n = strlen(s);
    char *out = (char*)malloc(n + 1);
    static const char map[256] = {
        ['A']='T', ['C']='G', ['G']='C', ['T']='A',
        ['a']='t', ['c']='g', ['g']='c', ['t']='a',
        ['N']='N', ['n']='n',
    };
    for (size_t i = 0; i < n; i++) {
        char c = map[(unsigned char)s[n - 1 - i]];
        out[i] = c ? c : 'N';
    }
    out[n] = '\0';
    return out;
}

static uint32_t g_rng = 0;
static void rng_seed(uint32_t s) { g_rng = s ? s : 1; }
static uint32_t rng_next(void) {
    uint32_t x = g_rng;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    g_rng = x;
    return x;
}
static void rand_seq(char *buf, int n) {
    static const char b[4] = {'A','C','G','T'};
    for (int i = 0; i < n; i++) buf[i] = b[rng_next() & 3];
    buf[n] = '\0';
}

#endif
