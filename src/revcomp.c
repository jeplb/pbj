#include "revcomp.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static char g_rc_table[256];
static int  g_rc_initialized = 0;

void pbj_revcomp_init(void) {
    if (g_rc_initialized) return;
    for (int i = 0; i < 256; i++) g_rc_table[i] = 'N';
    g_rc_table[(int)'A'] = 'T'; g_rc_table[(int)'T'] = 'A';
    g_rc_table[(int)'C'] = 'G'; g_rc_table[(int)'G'] = 'C';
    g_rc_table[(int)'a'] = 't'; g_rc_table[(int)'t'] = 'a';
    g_rc_table[(int)'c'] = 'g'; g_rc_table[(int)'g'] = 'c';
    g_rc_table[(int)'N'] = 'N'; g_rc_table[(int)'n'] = 'n';
    g_rc_table[(int)'U'] = 'A'; g_rc_table[(int)'u'] = 'a';
    g_rc_initialized = 1;
}

void pbj_revcomp(const pbj_read_t *src, pbj_read_t *dst) {
    uint32_t n = src->len;
    pbj_read_ensure_cap(dst, n + 1);
    dst->len = n;

    if (src->name) {
        size_t nlen = strlen(src->name);
        char *new_name = (char*)realloc(dst->name, nlen + 1);
        if (!new_name) {
            fprintf(stderr, "pbj: out of memory in pbj_revcomp name copy\n");
            abort();
        }
        dst->name = new_name;
        memcpy(dst->name, src->name, nlen + 1);
        dst->name_len = src->name_len;
    }

    const char *s = src->seq;
    const char *q = src->qual;
    char       *ds = dst->seq;
    char       *dq = dst->qual;
    for (uint32_t i = 0; i < n; i++) {
        ds[i] = g_rc_table[(unsigned char)s[n - 1 - i]];
        dq[i] = q[n - 1 - i];
    }
    ds[n] = '\0';
    dq[n] = '\0';
}
