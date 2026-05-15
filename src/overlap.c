#include "overlap.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* generation-counter scheme: instead of memset'ing the table to "empty"
   between pairs, bump cur_gen. a slot is occupied iff entry.gen ==
   cur_gen. on uint32 wrap we do a full reset. the gen field reuses the
   4 bytes of struct padding that would exist otherwise (uint64 + int32
   already round to 16). */
typedef struct {
    uint64_t kmer;
    int32_t  pos;
    uint32_t gen;
} kmer_entry_t;

typedef struct {
    kmer_entry_t *table;
    int           size;
    int           mask;
    uint32_t      cur_gen;
} kmer_index_t;

static inline int encode_base(int c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
    }
    return -1;
}

static inline uint32_t hash_u64(uint64_t x) {
    return (uint32_t)((x * 11400714819323198485ULL) >> 32);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) {
        fprintf(stderr, "pbj: out of memory (%zu bytes)\n", n);
        abort();
    }
    return p;
}

static void *xrealloc(void *q, size_t n) {
    void *p = realloc(q, n);
    if (!p) {
        fprintf(stderr, "pbj: out of memory (%zu bytes)\n", n);
        abort();
    }
    return p;
}

void pbj_overlap_init(pbj_workspace_t *ws, const pbj_params_t *p) {
    (void)p;
    kmer_index_t *ki = (kmer_index_t*)xmalloc(sizeof(*ki));
    ki->size  = 16384;
    ki->mask  = ki->size - 1;
    ki->table = (kmer_entry_t*)xmalloc((size_t)ki->size * sizeof(kmer_entry_t));
    /* gen 0 in every slot vs cur_gen starting at 1 means all slots
       initially read as empty without an O(size) zero pass. */
    for (int i = 0; i < ki->size; i++) ki->table[i].gen = 0;
    ki->cur_gen = 1;
    ws->kmer_index = ki;
    ws->candidates_cap = 256;
    ws->candidates = (int*)xmalloc((size_t)ws->candidates_cap * sizeof(int));
    ws->cand_bitset_bytes = 256;  /* enough for 2048 offsets, grows as needed */
    ws->cand_bitset = (uint8_t*)xmalloc((size_t)ws->cand_bitset_bytes);
}

void pbj_overlap_free(pbj_workspace_t *ws) {
    kmer_index_t *ki = (kmer_index_t*)ws->kmer_index;
    if (ki) {
        free(ki->table);
        free(ki);
        ws->kmer_index = NULL;
    }
    free(ws->candidates);
    ws->candidates = NULL;
    ws->candidates_cap = 0;
    free(ws->cand_bitset);
    ws->cand_bitset = NULL;
    ws->cand_bitset_bytes = 0;
}

int pbj_overlap_geometry(int s, int len_r1, int len_r2,
                          int min_overlap, int max_overlap,
                          int *r1_start, int *r2_start, int *olen) {
    int rs1, rs2, ol;
    if (s >= 0) {
        rs1 = s;
        rs2 = 0;
        int a = len_r1 - s;
        ol = (a < len_r2) ? a : len_r2;
    } else {
        rs1 = 0;
        rs2 = -s;
        int a = len_r2 + s;
        ol = (a < len_r1) ? a : len_r1;
    }
    if (ol < min_overlap) return 0;
    if (max_overlap > 0 && ol > max_overlap) return 0;
    *r1_start = rs1;
    *r2_start = rs2;
    *olen = ol;
    return 1;
}

int pbj_overlap_candidates(pbj_workspace_t *ws,
                            const pbj_read_t *r1,
                            const pbj_read_t *r2_rc,
                            const pbj_params_t *p) {
    kmer_index_t *ki = (kmer_index_t*)ws->kmer_index;
    int k  = p->kmer_size;
    int n1 = (int)r1->len;
    int n2 = (int)r2_rc->len;
    if (n1 < k || n2 < k) return 0;
    if (k <= 0 || k > 31) return 0;

    /* grow kmer table to keep load factor at most 0.25 */
    int needed = n1 + 1;
    int wanted = 16;
    while (wanted < 4 * needed) wanted <<= 1;
    if (ki->size < wanted) {
        free(ki->table);
        ki->table = (kmer_entry_t*)xmalloc((size_t)wanted * sizeof(kmer_entry_t));
        ki->size  = wanted;
        ki->mask  = wanted - 1;
        /* fresh allocation: zero gen so every slot starts as empty
           relative to any cur_gen >= 1. */
        for (int i = 0; i < ki->size; i++) ki->table[i].gen = 0;
    }
    /* mark every previous-batch entry as stale by bumping the generation.
       on uint32 wrap (extremely rare in practice; would require ~4
       billion pairs in one process) fall back to a full zero pass and
       restart at 1. */
    ki->cur_gen++;
    if (ki->cur_gen == 0) {
        for (int i = 0; i < ki->size; i++) ki->table[i].gen = 0;
        ki->cur_gen = 1;
    }

    uint64_t mask_k = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);

    /* index r1 kmers */
    int valid = 0;
    uint64_t kmer = 0;
    for (int i = 0; i < n1; i++) {
        int v = encode_base((unsigned char)r1->seq[i]);
        if (v < 0) { valid = 0; kmer = 0; continue; }
        kmer = ((kmer << 2) | (uint64_t)v) & mask_k;
        valid++;
        if (valid >= k) {
            int pos = i - k + 1;
            uint32_t h = hash_u64(kmer) & (uint32_t)ki->mask;
            while (ki->table[h].gen == ki->cur_gen) {
                h = (h + 1) & (uint32_t)ki->mask;
            }
            ki->table[h].kmer = kmer;
            ki->table[h].pos  = pos;
            ki->table[h].gen  = ki->cur_gen;
        }
    }

    int min_s = p->allow_dovetail ? -(n2 - p->min_overlap) : 0;
    int max_s = n1 - p->min_overlap;
    if (p->max_overlap > 0) {
        int lo = n1 - p->max_overlap;
        if (lo > min_s) min_s = lo;
    }
    if (min_s > max_s) return 0;

    /* dedup via a bitset over [min_s, max_s]. constant-time lookup vs the
       previous quadratic linear scan. */
    int span = max_s - min_s + 1;
    int bytes_needed = (span + 7) >> 3;
    if (ws->cand_bitset_bytes < bytes_needed) {
        ws->cand_bitset = (uint8_t*)xrealloc(ws->cand_bitset, (size_t)bytes_needed);
        ws->cand_bitset_bytes = bytes_needed;
    }
    memset(ws->cand_bitset, 0, (size_t)bytes_needed);

    int n_cand = 0;
    valid = 0; kmer = 0;
    for (int i = 0; i < n2; i++) {
        int v = encode_base((unsigned char)r2_rc->seq[i]);
        if (v < 0) { valid = 0; kmer = 0; continue; }
        kmer = ((kmer << 2) | (uint64_t)v) & mask_k;
        valid++;
        if (valid >= k) {
            int pos2 = i - k + 1;
            uint32_t h = hash_u64(kmer) & (uint32_t)ki->mask;
            while (ki->table[h].gen == ki->cur_gen) {
                if (ki->table[h].kmer == kmer) {
                    int s = ki->table[h].pos - pos2;
                    if (s >= min_s && s <= max_s) {
                        int bit_idx = s - min_s;
                        uint8_t mask = (uint8_t)(1u << (bit_idx & 7));
                        if (!(ws->cand_bitset[bit_idx >> 3] & mask)) {
                            ws->cand_bitset[bit_idx >> 3] |= mask;
                            if (n_cand >= ws->candidates_cap) {
                                ws->candidates_cap *= 2;
                                ws->candidates = (int*)xrealloc(ws->candidates,
                                    (size_t)ws->candidates_cap * sizeof(int));
                            }
                            ws->candidates[n_cand++] = s;
                        }
                    }
                }
                h = (h + 1) & (uint32_t)ki->mask;
            }
        }
    }
    return n_cand;
}
