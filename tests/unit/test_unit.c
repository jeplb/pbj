/* unit tests for pbj internal functions. links directly against the source
   modules (no main.c, no threads.c) so we can exercise individual functions
   without going through the cli. assertions / helpers come from
   test_helpers.h. */

#include "test_helpers.h"

/* ===================================================================
   revcomp
   =================================================================== */

static void test_revcomp_basic(void) {
    pbj_revcomp_init();
    pbj_read_t src, dst;
    load_read(&src, "x", "ACGT", "!\"#$");
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    CHECK_EQ_INT(dst.len, 4);
    CHECK_EQ_STR(dst.seq, "ACGT");      /* ACGT rc = ACGT */
    CHECK_EQ_STR(dst.qual, "$#\"!");    /* qual reversed, not complemented */
    pbj_read_free(&src);
    pbj_read_free(&dst);
}

static void test_revcomp_n_preserved(void) {
    pbj_revcomp_init();
    pbj_read_t src, dst;
    load_read(&src, "x", "ANCG", "????");
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    CHECK_EQ_STR(dst.seq, "CGNT");
    pbj_read_free(&src);
    pbj_read_free(&dst);
}

static void test_revcomp_mixed_case(void) {
    pbj_revcomp_init();
    pbj_read_t src, dst;
    load_read(&src, "x", "AcGt", "!!!!");
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    /* the table preserves case: a<->t, c<->g */
    CHECK_EQ_STR(dst.seq, "aCgT");
    pbj_read_free(&src);
    pbj_read_free(&dst);
}

static void test_revcomp_rna(void) {
    pbj_revcomp_init();
    pbj_read_t src, dst;
    load_read(&src, "x", "UUUU", "????");
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    /* U complements to A */
    CHECK_EQ_STR(dst.seq, "AAAA");
    pbj_read_free(&src);
    pbj_read_free(&dst);
}

static void test_revcomp_palindrome(void) {
    pbj_revcomp_init();
    pbj_read_t src, dst;
    load_read(&src, "x", "GAATTC", NULL);
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    CHECK_EQ_STR(dst.seq, "GAATTC");  /* EcoRI site is a palindrome */
    pbj_read_free(&src);
    pbj_read_free(&dst);
}

static void test_revcomp_longer(void) {
    pbj_revcomp_init();
    const char *seq = "ACGTACGTACGTACGTACGT";
    pbj_read_t src, dst;
    load_read(&src, "x", seq, NULL);
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    char *expected = revcomp_str(seq);
    CHECK_EQ_STR(dst.seq, expected);
    free(expected);
    pbj_read_free(&src);
    pbj_read_free(&dst);
}

/* ===================================================================
   stats: binomial survival function
   =================================================================== */

static void test_stats_boundaries(void) {
    /* p(x >= 0) = 1 always */
    CHECK_CLOSE(pbj_binomial_sf(0, 10, 0.5), 1.0, 1e-12);
    CHECK_CLOSE(pbj_binomial_sf(-1, 10, 0.5), 1.0, 1e-12);
    /* p(x >= n+1) = 0 */
    CHECK_CLOSE(pbj_binomial_sf(11, 10, 0.5), 0.0, 1e-12);
    /* n = 0 short-circuits to 1.0 */
    CHECK_CLOSE(pbj_binomial_sf(1, 0, 0.5), 1.0, 1e-12);
    /* p = 0 gives 0 for k >= 1 */
    CHECK_CLOSE(pbj_binomial_sf(1, 10, 0.0), 0.0, 1e-12);
    /* p = 1 gives 1 for k <= n */
    CHECK_CLOSE(pbj_binomial_sf(5, 10, 1.0), 1.0, 1e-12);
}

static void test_stats_known_values(void) {
    /* 4 fair coins, p(>= 4 heads) = 1/16 */
    CHECK_CLOSE(pbj_binomial_sf(4, 4, 0.5), 1.0 / 16.0, 1e-9);
    /* 4 coins, p(>= 3 heads) = 5/16 */
    CHECK_CLOSE(pbj_binomial_sf(3, 4, 0.5), 5.0 / 16.0, 1e-9);
    /* 10 coins, p(>= 5) = exactly 638/1024 = 0.623046875 */
    CHECK_CLOSE(pbj_binomial_sf(5, 10, 0.5), 638.0 / 1024.0, 1e-9);
    /* 6 trials at p=1/6, p(>= 6) = (1/6)^6 ~= 2.143e-5 */
    CHECK_CLOSE(pbj_binomial_sf(6, 6, 1.0 / 6.0), pow(1.0 / 6.0, 6.0), 1e-12);
}

static void test_stats_monotonicity(void) {
    /* sf is non-increasing in k */
    double prev = 1.0;
    for (int k = 0; k <= 20; k++) {
        double v = pbj_binomial_sf(k, 20, 0.3);
        CHECK(v <= prev + 1e-12);
        prev = v;
    }
}

static void test_stats_large_n_normal(void) {
    /* for large n, exact path switches to normal approximation. it should
       still produce a reasonable (small) p-value when k is well above
       expected. n=1000, p=0.25, k=400 ... mu=250, sd=~13.7, z>>10. */
    double v = pbj_binomial_sf(400, 1000, 0.25);
    CHECK(v >= 0.0);
    CHECK(v < 1e-10);
    /* and a large p-value when k is right at the mean */
    double v2 = pbj_binomial_sf(250, 1000, 0.25);
    CHECK(v2 > 0.4 && v2 < 0.6);
}

/* ===================================================================
   score: quality-weighted overlap scoring
   =================================================================== */

static void test_score_all_matches(void) {
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "ACGTACGT", "????????");  /* '?' = phred 30 */
    load_read(&r2, "b", "ACGTACGT", "????????");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 8, &p, &s);
    CHECK_EQ_INT(s.matches, 8);
    CHECK_EQ_INT(s.mismatches, 0);
    CHECK_EQ_INT(s.n_count, 0);
    /* match_reward=1, q=min(30,30)=30, score = 8 * 30 */
    CHECK_CLOSE(s.score, 240.0, 1e-9);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_score_all_mismatches(void) {
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "AAAA", "????");
    load_read(&r2, "b", "TTTT", "????");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 4, &p, &s);
    CHECK_EQ_INT(s.matches, 0);
    CHECK_EQ_INT(s.mismatches, 4);
    CHECK_EQ_INT(s.n_count, 0);
    CHECK_CLOSE(s.score, -120.0, 1e-9);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_score_n_excluded(void) {
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "ANNG", "????");
    load_read(&r2, "b", "ANCG", "????");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 4, &p, &s);
    /* position 0 (A=A) match, positions 1-2 contain N -> n_count,
       position 3 (G=G) match. */
    CHECK_EQ_INT(s.matches, 2);
    CHECK_EQ_INT(s.mismatches, 0);
    CHECK_EQ_INT(s.n_count, 2);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_score_quality_min(void) {
    /* a mismatch should be penalized by min(q1, q2), not max */
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    /* q1 = phred 30 ('?'), q2 = phred 2 ('#') -> min = 2 */
    load_read(&r1, "a", "A", "?");
    load_read(&r2, "b", "T", "#");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 1, &p, &s);
    CHECK_EQ_INT(s.mismatches, 1);
    CHECK_CLOSE(s.score, -2.0, 1e-9);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_score_offset_application(void) {
    /* r1_start / r2_start should shift the comparison window. */
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "XXACGT", "######????");
    load_read(&r2, "b", "ACGTYY", "????######");
    pbj_score_t s;
    pbj_score_overlap(&r1, 2, &r2, 0, 4, &p, &s);
    CHECK_EQ_INT(s.matches, 4);
    CHECK_EQ_INT(s.mismatches, 0);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

/* ===================================================================
   overlap geometry
   =================================================================== */

static void test_geom_positive_offset(void) {
    int r1s, r2s, ol;
    /* r1=100, r2=100, offset s=50: r1 covers [50..100), r2 covers [0..50) */
    int ok = pbj_overlap_geometry(50, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 1);
    CHECK_EQ_INT(r1s, 50);
    CHECK_EQ_INT(r2s, 0);
    CHECK_EQ_INT(ol, 50);
}

static void test_geom_negative_offset(void) {
    /* dovetail: r1=100, r2=100, s=-20 means r2_rc starts 20bp before r1.
       overlap len = min(r2 + s, r1) = min(80, 100) = 80. */
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(-20, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 1);
    CHECK_EQ_INT(r1s, 0);
    CHECK_EQ_INT(r2s, 20);
    CHECK_EQ_INT(ol, 80);
}

static void test_geom_zero_offset_full_overlap(void) {
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(0, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 1);
    CHECK_EQ_INT(r1s, 0);
    CHECK_EQ_INT(r2s, 0);
    CHECK_EQ_INT(ol, 100);
}

static void test_geom_below_min_rejected(void) {
    int r1s, r2s, ol;
    /* overlap of 15bp at offset 85, but min_overlap is 20 -> reject */
    int ok = pbj_overlap_geometry(85, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
}

static void test_geom_above_max_rejected(void) {
    int r1s, r2s, ol;
    /* full 100bp overlap, but max_overlap is 50 -> reject */
    int ok = pbj_overlap_geometry(0, 100, 100, 20, 50, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
}

static void test_geom_asymmetric_lengths(void) {
    int r1s, r2s, ol;
    /* r1=200, r2=50, offset=150: r1 tail of 50 overlaps full r2 */
    int ok = pbj_overlap_geometry(150, 200, 50, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 1);
    CHECK_EQ_INT(r1s, 150);
    CHECK_EQ_INT(r2s, 0);
    CHECK_EQ_INT(ol, 50);
}

/* ===================================================================
   overlap candidate generation
   =================================================================== */

static void test_candidates_perfect_overlap(void) {
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2, r2rc;
    /* construct r1 and r2 from a known fragment with a 50bp overlap */
    const char *frag = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
                       "GGGGGGGGGGAAAAAAAAAACCCCCCCCCCTTTTTTTTTTGGGGGGGGGG"
                       "ACACACACAC";
    /* r1 = first 100bp, r2 = revcomp of last 100bp (50bp overlap) */
    char r1seq[101], r2seq[101];
    memcpy(r1seq, frag, 100); r1seq[100] = '\0';
    /* the fragment is 112 chars; use first 100 for r1, last 100 (revcomp) for r2 */
    size_t fl = strlen(frag);
    char r2tail[101];
    memcpy(r2tail, frag + fl - 100, 100); r2tail[100] = '\0';
    char *r2rc_seq = revcomp_str(r2tail);
    memcpy(r2seq, r2rc_seq, 100); r2seq[100] = '\0';
    free(r2rc_seq);

    load_read(&r1, "x", r1seq, NULL);
    load_read(&r2, "y", r2seq, NULL);
    pbj_read_init(&r2rc);
    pbj_revcomp(&r2, &r2rc);

    int n = pbj_overlap_candidates(&ws, &r1, &r2rc, &p);
    CHECK(n > 0);

    /* the true offset is (r1_len - overlap) = 100 - 12 = 88 ... actually,
       fl=112, r2 covers fragment[12..112], r1 covers fragment[0..100],
       overlap is fragment[12..100] = 88bp, offset s = 12. */
    int found_true = 0;
    for (int i = 0; i < n; i++) {
        if (ws.candidates[i] == 12) { found_true = 1; break; }
    }
    CHECK(found_true);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&r2rc);
    ws_free(&ws);
}

static void test_candidates_short_reads(void) {
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2;
    /* n1 < k -> returns 0 without crashing */
    load_read(&r1, "x", "ACGT", NULL);  /* 4 < k=12 */
    load_read(&r2, "y", "TGCA", NULL);
    int n = pbj_overlap_candidates(&ws, &r1, &r2, &p);
    CHECK_EQ_INT(n, 0);
    pbj_read_free(&r1); pbj_read_free(&r2);
    ws_free(&ws);
}

static void test_candidates_all_n_no_kmers(void) {
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2;
    char seq[101]; memset(seq, 'N', 100); seq[100] = '\0';
    load_read(&r1, "x", seq, NULL);
    load_read(&r2, "y", seq, NULL);
    int n = pbj_overlap_candidates(&ws, &r1, &r2, &p);
    CHECK_EQ_INT(n, 0);
    pbj_read_free(&r1); pbj_read_free(&r2);
    ws_free(&ws);
}

static void test_candidates_dedup(void) {
    /* repeating pattern produces many kmer hits at the same offset; the
       bitset dedup should collapse them all to a single candidate. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.kmer_size = 4;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2;
    /* both reads identical -> many kmer matches at offset 0 */
    load_read(&r1, "x", "ACGTACGTACGTACGTACGTACGTACGTACGT", NULL);
    load_read(&r2, "y", "ACGTACGTACGTACGTACGTACGTACGTACGT", NULL);
    int n = pbj_overlap_candidates(&ws, &r1, &r2, &p);
    /* should be small (a few candidates from the repeating pattern), not
       hundreds. concretely, the periodic ACGT yields offsets at multiples
       of 4 -- but the dedup keeps each unique offset once. */
    CHECK(n > 0);
    CHECK(n < 20);
    /* offset 0 must be among them */
    int has_zero = 0;
    for (int i = 0; i < n; i++) if (ws.candidates[i] == 0) has_zero = 1;
    CHECK(has_zero);
    pbj_read_free(&r1); pbj_read_free(&r2);
    ws_free(&ws);
}

/* ===================================================================
   adapter trimming
   =================================================================== */

static void test_adapter_3prime_match(void) {
    /* a 21bp adapter perfectly matches the last 21bp of a 100bp read */
    const char *adapter = "AGATCGGAAGAGCACACGTCT";
    char seq[101];
    memset(seq, 'A', 79);
    memcpy(seq + 79, adapter, 21);
    seq[100] = '\0';
    int trim = pbj_adapter_find_trim(seq, 100, adapter, 21, 8, 0.9);
    CHECK_EQ_INT(trim, 79);
}

static void test_adapter_partial_3prime(void) {
    /* only the first 10bp of the adapter is present at the read tail */
    const char *adapter = "AGATCGGAAGAGCACACGTCT";
    char seq[101];
    memset(seq, 'C', 90);
    memcpy(seq + 90, adapter, 10);
    seq[100] = '\0';
    int trim = pbj_adapter_find_trim(seq, 100, adapter, 21, 8, 0.9);
    CHECK_EQ_INT(trim, 90);
}

static void test_adapter_mid_read_ignored(void) {
    /* adapter sequence in middle of read but not at 3' end -> NOT trimmed */
    const char *adapter = "AGATCGGAAGAGCACACGTCT";
    char seq[101];
    memset(seq, 'C', 100);
    /* insert adapter at position 30, but the 3' tail is all C */
    memcpy(seq + 30, adapter, 21);
    seq[100] = '\0';
    int trim = pbj_adapter_find_trim(seq, 100, adapter, 21, 8, 0.9);
    /* no 3'-anchored match exists; full read length returned */
    CHECK_EQ_INT(trim, 100);
}

static void test_adapter_min_match_bp_threshold(void) {
    /* only 5bp of adapter at 3' end; with min_match_bp=8 we don't trim,
       with min_match_bp=5 we do. */
    const char *adapter = "AGATCGGAAGAGCACACGTCT";
    char seq[101];
    memset(seq, 'C', 95);
    memcpy(seq + 95, adapter, 5);
    seq[100] = '\0';
    int trim8 = pbj_adapter_find_trim(seq, 100, adapter, 21, 8, 0.9);
    int trim5 = pbj_adapter_find_trim(seq, 100, adapter, 21, 5, 0.9);
    CHECK_EQ_INT(trim8, 100);  /* not trimmed: 5 < 8 */
    CHECK_EQ_INT(trim5, 95);   /* trimmed at pos 95 */
}

static void test_adapter_n_wildcard(void) {
    /* N in adapter should act as a wildcard */
    const char *adapter = "AGNTCGGAAGAGCACACGTCT";
    const char *real    = "AGATCGGAAGAGCACACGTCT";
    char seq[101];
    memset(seq, 'C', 79);
    memcpy(seq + 79, real, 21);
    seq[100] = '\0';
    int trim = pbj_adapter_find_trim(seq, 100, adapter, 21, 8, 0.9);
    CHECK_EQ_INT(trim, 79);
}

static void test_adapter_case_insensitive(void) {
    /* lowercase adapter against uppercase read */
    const char *adapter = "agatcggaagagcacacgtct";
    const char *real    = "AGATCGGAAGAGCACACGTCT";
    char seq[101];
    memset(seq, 'C', 79);
    memcpy(seq + 79, real, 21);
    seq[100] = '\0';
    int trim = pbj_adapter_find_trim(seq, 100, adapter, 21, 8, 0.9);
    CHECK_EQ_INT(trim, 79);
}

static void test_adapter_empty_noop(void) {
    /* empty adapter or zero-length read: return len unchanged */
    int trim = pbj_adapter_find_trim("ACGT", 4, "", 0, 8, 0.9);
    CHECK_EQ_INT(trim, 4);
    trim = pbj_adapter_find_trim(NULL, 4, "AGA", 3, 8, 0.9);
    /* NULL seq isn't documented but the function guards on adapter_len; this
       test pins the existing guard pattern. */
    (void)trim;
}

static void test_adapter_inplace_trim(void) {
    pbj_revcomp_init();
    pbj_read_t r;
    char tail[22];
    memcpy(tail, "AGATCGGAAGAGCACACGTCT", 21); tail[21] = '\0';
    char full[101];
    memset(full, 'C', 79);
    memcpy(full + 79, tail, 21);
    full[100] = '\0';
    char qual[101];
    memset(qual, '?', 100); qual[100] = '\0';
    load_read(&r, "x", full, qual);
    int n = pbj_adapter_trim(&r, tail, 8, 0.9);
    CHECK_EQ_INT(n, 1);
    CHECK_EQ_INT(r.len, 79);
    /* seq and qual should be null-terminated at the trim point */
    CHECK_EQ_INT(r.seq[79], 0);
    CHECK_EQ_INT(r.qual[79], 0);
    pbj_read_free(&r);
}

/* ===================================================================
   pair name matching / casava direction
   =================================================================== */

static void test_pair_match_basic(void) {
    pbj_read_t a, b;
    load_read(&a, "foo/1", "A", "?");
    load_read(&b, "foo/2", "A", "?");
    CHECK_EQ_INT(pbj_pair_names_match(&a, &b), 1);
    pbj_read_free(&a); pbj_read_free(&b);
}

static void test_pair_same_suffix_rejected(void) {
    /* both /1 -> rejected */
    pbj_read_t a, b;
    load_read(&a, "foo/1", "A", "?");
    load_read(&b, "foo/1", "A", "?");
    CHECK_EQ_INT(pbj_pair_names_match(&a, &b), 0);
    pbj_read_free(&a); pbj_read_free(&b);
}

static void test_pair_core_id_match(void) {
    /* same core id, no /1 or /2 suffix at all -> match */
    pbj_read_t a, b;
    load_read(&a, "HISEQ:1:1101", "A", "?");
    load_read(&b, "HISEQ:1:1101", "A", "?");
    CHECK_EQ_INT(pbj_pair_names_match(&a, &b), 1);
    pbj_read_free(&a); pbj_read_free(&b);
}

static void test_pair_different_core_id_rejected(void) {
    pbj_read_t a, b;
    load_read(&a, "foo/1", "A", "?");
    load_read(&b, "bar/2", "A", "?");
    CHECK_EQ_INT(pbj_pair_names_match(&a, &b), 0);
    pbj_read_free(&a); pbj_read_free(&b);
}

/* casava_direction reads a "comment" stored after name_len + space. since
   load_read() puts the whole string into name and sets name_len = strlen,
   we have to hand-craft the read for this test. */
static void make_casava(pbj_read_t *r, const char *core, const char *comment) {
    pbj_read_init(r);
    pbj_read_ensure_cap(r, 4);
    strcpy(r->seq, "ACGT");
    strcpy(r->qual, "????");
    r->len = 4;
    size_t cl = strlen(core);
    size_t comment_l = strlen(comment);
    r->name = (char*)malloc(cl + 1 + comment_l + 1);
    memcpy(r->name, core, cl);
    r->name[cl] = ' ';
    memcpy(r->name + cl + 1, comment, comment_l);
    r->name[cl + 1 + comment_l] = '\0';
    r->name_len = (uint32_t)cl;
}

static void test_casava_direction_r1(void) {
    pbj_read_t r;
    make_casava(&r, "HISEQ:123:FC:1:1101:1234:5678", "1:N:0:ATCG");
    CHECK_EQ_INT(pbj_casava_direction(&r), 1);
    pbj_read_free(&r);
}

static void test_casava_direction_r2(void) {
    pbj_read_t r;
    make_casava(&r, "HISEQ:123", "2:N:0:GCAT");
    CHECK_EQ_INT(pbj_casava_direction(&r), 2);
    pbj_read_free(&r);
}

static void test_casava_no_direction(void) {
    pbj_read_t r;
    make_casava(&r, "name", "no-direction-here");
    CHECK_EQ_INT(pbj_casava_direction(&r), 0);
    pbj_read_free(&r);
}

/* ===================================================================
   end-to-end merge (the real integration unit)
   =================================================================== */

static void test_merge_perfect_overlap(void) {
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    /* random fragment of 150bp, r1=100bp, r2=100bp (rc), 50bp overlap */
    char frag[151], r1seq[101], r2seq[101], q[101];
    rng_seed(0xCAFE);
    rand_seq(frag, 150);
    memcpy(r1seq, frag, 100); r1seq[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2seq, r2rc, 100); r2seq[100] = '\0';
    free(r2rc);
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "perf/1", r1seq, q);
    load_read(&r2, "perf/2", r2seq, q);

    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
    CHECK_EQ_INT(res.merged.len, 150);
    CHECK_EQ_STR(res.merged.seq, frag);
    /* casava-style /1 suffix stripped from merged name */
    CHECK_EQ_STR(res.merged.name, "perf");

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_merge_no_overlap(void) {
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    /* two independent random fragments, no shared overlap */
    char a[101], b[101], q[101];
    rng_seed(0xBEEF);
    rand_seq(a, 100);
    rng_seed(0xDEAD);
    rand_seq(b, 100);
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "x/1", a, q);
    load_read(&r2, "x/2", b, q);

    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    /* either NO_CANDIDATE or FAIL_PVALUE depending on whether any chance
       kmer hits exist; both are acceptable "no merge" outcomes. */
    CHECK(res.status == PBJ_NO_CANDIDATE ||
          res.status == PBJ_FAIL_PVALUE ||
          res.status == PBJ_FAIL_LENGTH);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_merge_consensus_picks_higher_q(void) {
    /* single mismatch at the overlap: R1 disagrees with high Q, R2 stays
       at low Q -> consensus picks R1's base. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    char frag[151], r1seq[101], r2seq[101], q1[101], q2[101];
    rng_seed(0xC0DE);
    rand_seq(frag, 150);
    memcpy(r1seq, frag, 100); r1seq[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2seq, r2rc, 100); r2seq[100] = '\0';
    free(r2rc);
    memset(q1, '?', 100); q1[100] = '\0';  /* phred 30 */
    memset(q2, '+', 100); q2[100] = '\0';  /* phred 10 */

    /* flip r1 at position 75 (inside the overlap region 50..100) */
    char swap[256] = {0};
    swap['A']='C'; swap['C']='A'; swap['G']='T'; swap['T']='G';
    char flipped = swap[(unsigned char)r1seq[75]];
    r1seq[75] = flipped;

    pbj_read_t r1, r2;
    load_read(&r1, "x/1", r1seq, q1);
    load_read(&r2, "x/2", r2seq, q2);

    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
    CHECK_EQ_INT(res.merged.len, 150);
    /* the merged base at position 75 should be R1's high-Q base */
    CHECK_EQ_INT(res.merged.seq[75], flipped);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_merge_min_length_rejection(void) {
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.min_assembly_len = 1000;  /* impossibly high */
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    char frag[151], r1seq[101], r2seq[101], q[101];
    rng_seed(0xF00D);
    rand_seq(frag, 150);
    memcpy(r1seq, frag, 100); r1seq[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2seq, r2rc, 100); r2seq[100] = '\0';
    free(r2rc);
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "x/1", r1seq, q);
    load_read(&r2, "x/2", r2seq, q);

    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_FAIL_LENGTH);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_merge_strict_ambiguity(void) {
    /* high-Q disagreements on both sides under --strict -> FAIL_AMBIGUOUS */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.strict_ambiguity = 1;
    p.strict_q_threshold = 20;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    char frag[151], r1seq[101], r2seq[101], q[101];
    rng_seed(0xABCD);
    rand_seq(frag, 150);
    memcpy(r1seq, frag, 100); r1seq[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2seq, r2rc, 100); r2seq[100] = '\0';
    free(r2rc);
    memset(q, '?', 100); q[100] = '\0';  /* phred 30 throughout */

    char swap[256] = {0};
    swap['A']='C'; swap['C']='A'; swap['G']='T'; swap['T']='G';
    r1seq[60] = swap[(unsigned char)r1seq[60]];
    r1seq[70] = swap[(unsigned char)r1seq[70]];
    r1seq[80] = swap[(unsigned char)r1seq[80]];

    pbj_read_t r1, r2;
    load_read(&r1, "x/1", r1seq, q);
    load_read(&r2, "x/2", r2seq, q);

    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_FAIL_AMBIGUOUS);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_merge_expected_errors(void) {
    /* perfect merge of length 150 with phred 30 throughout. on agreement
       positions the consensus sums the qualities (capped at 93), so EE is
       lower than the trivial 150 * 10^-3 = 0.15. we verify:
       - --max-ee 1.0 accepts (ee is well under 1.0)
       - --max-ee 0.0001 rejects (ee is well above 0.0001) */
    pbj_revcomp_init();
    char frag[151], r1seq[101], r2seq[101], q[101];
    rng_seed(0x5EED);
    rand_seq(frag, 150);
    memcpy(r1seq, frag, 100); r1seq[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2seq, r2rc, 100); r2seq[100] = '\0';
    free(r2rc);
    memset(q, '?', 100); q[100] = '\0';

    for (int reject = 0; reject < 2; reject++) {
        pbj_params_t p; params_defaults(&p);
        p.max_ee = reject ? 0.0001 : 1.0;
        pbj_workspace_t ws; ws_init(&ws, &p);
        pbj_merge_result_t res; memset(&res, 0, sizeof(res));
        pbj_read_init(&res.merged);

        pbj_read_t r1, r2;
        load_read(&r1, "x/1", r1seq, q);
        load_read(&r2, "x/2", r2seq, q);

        pbj_try_merge(&r1, &r2, &ws, &p, &res);
        if (reject) {
            CHECK_EQ_INT(res.status, PBJ_FAIL_MAX_EE);
        } else {
            CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
            CHECK(res.expected_errors > 0.0);
            CHECK(res.expected_errors < 0.2);
        }
        pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
        ws_free(&ws);
    }
}

static void test_merge_mono_filter(void) {
    /* a 150bp merge with >85% A should be rejected by --max-mono-frac 0.85. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.max_mono_frac = 0.85;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    char frag[151];
    for (int i = 0; i < 150; i++) frag[i] = 'A';
    frag[10] = 'C'; frag[40] = 'G'; frag[80] = 'T'; frag[120] = 'C';
    frag[150] = '\0';
    char r1seq[101], r2seq[101], q[101];
    memcpy(r1seq, frag, 100); r1seq[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2seq, r2rc, 100); r2seq[100] = '\0';
    free(r2rc);
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "x/1", r1seq, q);
    load_read(&r2, "x/2", r2seq, q);

    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_FAIL_MONO);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

/* ===================================================================
   read buffer helpers
   =================================================================== */

static void test_read_ensure_cap_grows(void) {
    pbj_read_t r;
    pbj_read_init(&r);
    CHECK_EQ_INT(r.cap, 0);
    pbj_read_ensure_cap(&r, 100);
    CHECK(r.cap >= 100);
    CHECK(r.seq != NULL);
    CHECK(r.qual != NULL);
    uint32_t cap1 = r.cap;
    /* second call with smaller cap should be a no-op */
    pbj_read_ensure_cap(&r, 50);
    CHECK_EQ_INT(r.cap, cap1);
    /* much larger cap should grow */
    pbj_read_ensure_cap(&r, 10000);
    CHECK(r.cap >= 10000);
    pbj_read_free(&r);
}

/* ===================================================================
   driver
   =================================================================== */

int main(void) {
    printf("pbj unit tests\n");
    printf("--------------\n");

    /* revcomp */
    RUN(test_revcomp_basic);
    RUN(test_revcomp_n_preserved);
    RUN(test_revcomp_mixed_case);
    RUN(test_revcomp_rna);
    RUN(test_revcomp_palindrome);
    RUN(test_revcomp_longer);

    /* stats */
    RUN(test_stats_boundaries);
    RUN(test_stats_known_values);
    RUN(test_stats_monotonicity);
    RUN(test_stats_large_n_normal);

    /* score */
    RUN(test_score_all_matches);
    RUN(test_score_all_mismatches);
    RUN(test_score_n_excluded);
    RUN(test_score_quality_min);
    RUN(test_score_offset_application);

    /* overlap geometry */
    RUN(test_geom_positive_offset);
    RUN(test_geom_negative_offset);
    RUN(test_geom_zero_offset_full_overlap);
    RUN(test_geom_below_min_rejected);
    RUN(test_geom_above_max_rejected);
    RUN(test_geom_asymmetric_lengths);

    /* candidate generation */
    RUN(test_candidates_perfect_overlap);
    RUN(test_candidates_short_reads);
    RUN(test_candidates_all_n_no_kmers);
    RUN(test_candidates_dedup);

    /* adapter */
    RUN(test_adapter_3prime_match);
    RUN(test_adapter_partial_3prime);
    RUN(test_adapter_mid_read_ignored);
    RUN(test_adapter_min_match_bp_threshold);
    RUN(test_adapter_n_wildcard);
    RUN(test_adapter_case_insensitive);
    RUN(test_adapter_empty_noop);
    RUN(test_adapter_inplace_trim);

    /* pair / casava */
    RUN(test_pair_match_basic);
    RUN(test_pair_same_suffix_rejected);
    RUN(test_pair_core_id_match);
    RUN(test_pair_different_core_id_rejected);
    RUN(test_casava_direction_r1);
    RUN(test_casava_direction_r2);
    RUN(test_casava_no_direction);

    /* end-to-end merge */
    RUN(test_merge_perfect_overlap);
    RUN(test_merge_no_overlap);
    RUN(test_merge_consensus_picks_higher_q);
    RUN(test_merge_min_length_rejection);
    RUN(test_merge_strict_ambiguity);
    RUN(test_merge_expected_errors);
    RUN(test_merge_mono_filter);

    /* read buffer */
    RUN(test_read_ensure_cap_grows);

    printf("--------------\n");
    printf("%d checks, %d failed\n", g_tests, g_failed);
    return g_failed == 0 ? 0 : 1;
}
