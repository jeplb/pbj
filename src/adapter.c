#include "adapter.h"
#include <string.h>

static inline char up(char c) {
    return (c >= 'a' && c <= 'z') ? (char)(c - 32) : c;
}

/* find the leftmost 3'-anchored position in seq where a partial-prefix
   match of adapter reaches the requested match fraction. "3'-anchored"
   means the candidate adapter occurrence must extend all the way to the
   read end (pos + eff == len), i.e., we only consider adapter readthrough
   at the read tail. this is the standard semantic in cutadapt / fastp /
   bbduk for a 3' adapter and rules out random mid-read prefix collisions.
   the candidate range is therefore [max(0, len - adapter_len),
   len - min_match_bp]; we walk left to right and pick the first hit,
   which gives the longest acceptable match (closest to the full adapter)
   and trims the most contamination off the read. */
int pbj_adapter_find_trim(const char *seq, int len,
                          const char *adapter, int adapter_len,
                          int min_match_bp, double min_match_frac) {
    if (!adapter || adapter_len <= 0 || len <= 0) return len;
    if (min_match_bp <= 0)                          min_match_bp = 5;
    if (min_match_frac <= 0.0 || min_match_frac > 1.0) min_match_frac = 0.9;
    if (min_match_bp > len) return len;

    int start = (adapter_len < len) ? (len - adapter_len) : 0;
    int end_excl = len - min_match_bp + 1;
    if (end_excl > len) end_excl = len;

    for (int pos = start; pos < end_excl; pos++) {
        int eff = len - pos;
        int matches = 0;
        for (int i = 0; i < eff; i++) {
            /* reads are uppercased at parse time, but adapter strings come
               from argv unmodified. uppercase per-byte at compare time so
               -a agatcgg and -a AGATCGG behave identically. */
            char b1 = up(seq[pos + i]);
            char b2 = up(adapter[i]);
            if (b1 == 'N' || b2 == 'N' || b1 == b2) matches++;
        }
        if ((double)matches >= (double)eff * min_match_frac) {
            return pos;
        }
    }
    return len;
}

int pbj_adapter_trim(pbj_read_t *r,
                     const char *adapter,
                     int min_match_bp, double min_match_frac) {
    if (!adapter || !r->seq) return 0;
    int alen = (int)strlen(adapter);
    int new_len = pbj_adapter_find_trim(r->seq, (int)r->len, adapter, alen,
                                         min_match_bp, min_match_frac);
    if (new_len < (int)r->len) {
        r->len = (uint32_t)new_len;
        r->seq[new_len]  = '\0';
        r->qual[new_len] = '\0';
        return 1;
    }
    return 0;
}
