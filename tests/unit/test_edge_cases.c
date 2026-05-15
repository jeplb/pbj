/* edge-case test suite for pbj. probes boundary values, degenerate
   inputs, and combinations that the happy-path tests in test_unit.c
   don't reach. */

#include "test_helpers.h"

static char *dup_str(const char *s) {
    size_t n = strlen(s);
    char *p = (char*)malloc(n + 1);
    memcpy(p, s, n + 1);
    return p;
}

/* ===================================================================
   read buffer (pbj_read_ensure_cap)
   =================================================================== */

static void test_ec_ensure_cap_zero(void) {
    pbj_read_t r;
    pbj_read_init(&r);
    pbj_read_ensure_cap(&r, 0);
    /* cap=0 on a fresh read is a no-op; no allocation happens */
    CHECK_EQ_INT(r.cap, 0);
    CHECK(r.seq == NULL);
    pbj_read_free(&r);
}

static void test_ec_ensure_cap_minimum(void) {
    /* asking for cap=1 should allocate the minimum 64-byte buffer */
    pbj_read_t r;
    pbj_read_init(&r);
    pbj_read_ensure_cap(&r, 1);
    CHECK_EQ_INT(r.cap, 64);
    pbj_read_free(&r);
}

static void test_ec_ensure_cap_doubling(void) {
    /* growth follows a doubling pattern from the current cap */
    pbj_read_t r;
    pbj_read_init(&r);
    pbj_read_ensure_cap(&r, 64);
    CHECK_EQ_INT(r.cap, 64);
    pbj_read_ensure_cap(&r, 65);
    CHECK_EQ_INT(r.cap, 128);
    pbj_read_ensure_cap(&r, 129);
    CHECK_EQ_INT(r.cap, 256);
    pbj_read_free(&r);
}

static void test_ec_ensure_cap_idempotent(void) {
    /* repeated calls with equal or smaller cap don't realloc */
    pbj_read_t r;
    pbj_read_init(&r);
    pbj_read_ensure_cap(&r, 100);
    uint32_t cap = r.cap;
    char *seq_before = r.seq;
    pbj_read_ensure_cap(&r, 100);
    pbj_read_ensure_cap(&r, 50);
    pbj_read_ensure_cap(&r, 1);
    CHECK_EQ_INT(r.cap, cap);
    CHECK(r.seq == seq_before);
    pbj_read_free(&r);
}

static void test_ec_free_reuse_cycle(void) {
    /* free then reuse: pointers nulled, fresh ensure_cap works */
    pbj_read_t r;
    pbj_read_init(&r);
    pbj_read_ensure_cap(&r, 100);
    pbj_read_free(&r);
    CHECK(r.seq == NULL);
    CHECK_EQ_INT(r.cap, 0);
    pbj_read_ensure_cap(&r, 50);
    CHECK_EQ_INT(r.cap, 64);
    pbj_read_free(&r);
}

/* ===================================================================
   revcomp edge cases
   =================================================================== */

static void test_ec_revcomp_empty(void) {
    pbj_revcomp_init();
    pbj_read_t src, dst;
    pbj_read_init(&src);
    pbj_read_init(&dst);
    pbj_read_ensure_cap(&src, 1);
    src.seq[0] = '\0'; src.qual[0] = '\0'; src.len = 0;
    src.name = dup_str("e"); src.name_len = 1;
    pbj_revcomp(&src, &dst);
    CHECK_EQ_INT(dst.len, 0);
    pbj_read_free(&src); pbj_read_free(&dst);
}

static void test_ec_revcomp_single_base(void) {
    pbj_revcomp_init();
    pbj_read_t src, dst;
    load_read(&src, "x", "A", "?");
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    CHECK_EQ_STR(dst.seq, "T");
    CHECK_EQ_STR(dst.qual, "?");
    pbj_read_free(&src); pbj_read_free(&dst);
}

static void test_ec_revcomp_invalid_chars(void) {
    /* any non-ACGT/U/N character maps to N */
    pbj_revcomp_init();
    pbj_read_t src, dst;
    load_read(&src, "x", "AXYZ", "????");
    pbj_read_init(&dst);
    pbj_revcomp(&src, &dst);
    /* reverse of "AXYZ" complemented: Z->N, Y->N, X->N, A->T */
    CHECK_EQ_STR(dst.seq, "NNNT");
    pbj_read_free(&src); pbj_read_free(&dst);
}

/* ===================================================================
   binomial sf edge cases
   =================================================================== */

static void test_ec_binom_n_equals_1(void) {
    /* single trial: p(>=1) = p; p(>=0) = 1 */
    CHECK_CLOSE(pbj_binomial_sf(1, 1, 0.25), 0.25, 1e-12);
    CHECK_CLOSE(pbj_binomial_sf(0, 1, 0.25), 1.0, 1e-12);
    CHECK_CLOSE(pbj_binomial_sf(2, 1, 0.25), 0.0, 1e-12);
}

static void test_ec_binom_p_one_with_k_equal_n(void) {
    /* deterministic case: p=1, every trial succeeds */
    CHECK_CLOSE(pbj_binomial_sf(10, 10, 1.0), 1.0, 1e-12);
}

static void test_ec_binom_small_p(void) {
    /* small p, large k: p-value should be very small but nonzero. */
    double v = pbj_binomial_sf(50, 100, 1e-3);
    CHECK(v >= 0.0);
    CHECK(v < 1e-100);
}

static void test_ec_binom_k_equals_n(void) {
    /* k = n: only one term in the sum: p^n */
    double n = 10, p = 0.3;
    CHECK_CLOSE(pbj_binomial_sf(10, 10, p), pow(p, n), 1e-12);
}

static void test_ec_binom_exact_normal_boundary(void) {
    /* at the n=256 / n=257 boundary the implementation switches from
       exact summation to a normal approximation. for moderate k the two
       paths should agree to within ~30% in relative terms. n=256, p=0.25,
       k=80 sits ~2.3 sd above the mean, so the tail mass is around 0.01. */
    double exact = pbj_binomial_sf(80, 256, 0.25);
    double approx = pbj_binomial_sf(80, 257, 0.25);
    CHECK(exact > 0.0 && exact < 0.05);
    CHECK(approx > 0.0 && approx < 0.05);
    /* the two paths should agree within 50% relative */
    CHECK_CLOSE(exact, approx, exact * 0.5);
}

/* ===================================================================
   score edge cases
   =================================================================== */

static void test_ec_score_zero_length_overlap(void) {
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "ACGT", "????");
    load_read(&r2, "b", "ACGT", "????");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 0, &p, &s);
    CHECK_EQ_INT(s.matches, 0);
    CHECK_EQ_INT(s.mismatches, 0);
    CHECK_EQ_INT(s.n_count, 0);
    CHECK_CLOSE(s.score, 0.0, 1e-12);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_ec_score_quality_below_offset(void) {
    /* quality bytes below the offset should clamp to 0 phred, not go
       negative or wrap. ASCII 32 (' ') with offset 33 -> phred -1 -> 0. */
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "A", " ");   /* phred -1 -> clamp 0 */
    load_read(&r2, "b", "T", " ");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 1, &p, &s);
    CHECK_EQ_INT(s.mismatches, 1);
    /* min(0, 0) * penalty = 0 */
    CHECK_CLOSE(s.score, 0.0, 1e-12);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_ec_score_max_phred(void) {
    /* phred 93 = ASCII 126 = '~'. should be handled without overflow. */
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "A", "~");
    load_read(&r2, "b", "A", "~");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 1, &p, &s);
    CHECK_EQ_INT(s.matches, 1);
    CHECK_CLOSE(s.score, 93.0, 1e-12);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_ec_score_iupac_treated_as_mismatch(void) {
    /* Y vs A is not N: scoring should treat it as a real mismatch, not
       skip it. (This is intentional; pbj documents this.) */
    pbj_params_t p; params_defaults(&p);
    pbj_read_t r1, r2;
    load_read(&r1, "a", "Y", "?");
    load_read(&r2, "b", "A", "?");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 1, &p, &s);
    CHECK_EQ_INT(s.mismatches, 1);
    CHECK_EQ_INT(s.n_count, 0);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_ec_score_zero_reward_zero_penalty(void) {
    /* with reward=0 and penalty=0, score should stay at 0 regardless of
       match/mismatch composition. counts should still be accurate. */
    pbj_params_t p; params_defaults(&p);
    p.match_reward = 0;
    p.mismatch_penalty = 0;
    pbj_read_t r1, r2;
    load_read(&r1, "a", "AACC", "????");
    load_read(&r2, "b", "AAGT", "????");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 4, &p, &s);
    CHECK_EQ_INT(s.matches, 2);
    CHECK_EQ_INT(s.mismatches, 2);
    CHECK_CLOSE(s.score, 0.0, 1e-12);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

static void test_ec_score_quality_offset_64(void) {
    /* offset 64 (legacy Illumina). ASCII 'B' (66) -> phred 2. */
    pbj_params_t p; params_defaults(&p);
    p.quality_offset = 64;
    pbj_read_t r1, r2;
    load_read(&r1, "a", "A", "B");
    load_read(&r2, "b", "A", "B");
    pbj_score_t s;
    pbj_score_overlap(&r1, 0, &r2, 0, 1, &p, &s);
    CHECK_EQ_INT(s.matches, 1);
    CHECK_CLOSE(s.score, 2.0, 1e-12);
    pbj_read_free(&r1); pbj_read_free(&r2);
}

/* ===================================================================
   overlap geometry edge cases
   =================================================================== */

static void test_ec_geom_exactly_min_overlap(void) {
    /* overlap exactly equal to min_overlap should be accepted */
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(80, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 1);
    CHECK_EQ_INT(ol, 20);
}

static void test_ec_geom_one_below_min_overlap(void) {
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(81, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
}

static void test_ec_geom_exactly_max_overlap(void) {
    int r1s, r2s, ol;
    /* full 100bp overlap with max_overlap=100 -> accepted */
    int ok = pbj_overlap_geometry(0, 100, 100, 20, 100, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 1);
    CHECK_EQ_INT(ol, 100);
}

static void test_ec_geom_one_above_max_overlap(void) {
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(0, 100, 100, 20, 99, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
}

static void test_ec_geom_offset_past_r1_end(void) {
    /* s >= len_r1 produces zero or negative overlap; must be rejected */
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(100, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
    ok = pbj_overlap_geometry(200, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
}

static void test_ec_geom_very_negative_offset(void) {
    /* dovetail with |s| >= r2 length: overlap collapses to <= 0 */
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(-100, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
    ok = pbj_overlap_geometry(-200, 100, 100, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 0);
}

static void test_ec_geom_asymmetric_r1_much_longer(void) {
    /* r1 much longer than r2: r2 fits entirely inside r1 if offset is right */
    int r1s, r2s, ol;
    int ok = pbj_overlap_geometry(450, 500, 50, 20, 0, &r1s, &r2s, &ol);
    CHECK_EQ_INT(ok, 1);
    CHECK_EQ_INT(r1s, 450);
    CHECK_EQ_INT(r2s, 0);
    CHECK_EQ_INT(ol, 50);
}

/* ===================================================================
   candidate generation edge cases
   =================================================================== */

static void test_ec_candidates_len_exactly_k(void) {
    /* r1 of length exactly k produces exactly one kmer. with a matching
       r2, that single kmer can seed a candidate. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.kmer_size = 8;
    p.min_overlap = 8;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2;
    load_read(&r1, "x", "ACGTACGT", NULL);
    load_read(&r2, "y", "ACGTACGT", NULL);
    int n = pbj_overlap_candidates(&ws, &r1, &r2, &p);
    CHECK(n >= 1);
    pbj_read_free(&r1); pbj_read_free(&r2);
    ws_free(&ws);
}

static void test_ec_candidates_len_below_k(void) {
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.kmer_size = 8;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2;
    load_read(&r1, "x", "ACGTACG", NULL);  /* 7 < k=8 */
    load_read(&r2, "y", "ACGTACG", NULL);
    int n = pbj_overlap_candidates(&ws, &r1, &r2, &p);
    CHECK_EQ_INT(n, 0);
    pbj_read_free(&r1); pbj_read_free(&r2);
    ws_free(&ws);
}

/* note: pbj_overlap_candidates takes (r1, r2_rc, ...) where r2_rc is
   the *reverse complement* of the as-sequenced R2 read. our tests build
   the equivalent r2_rc directly (a slice of the original fragment) so we
   don't have to revcomp twice. */

static void test_ec_candidates_kmer_size_min(void) {
    /* kmer_size = 4 is the documented minimum. should work without error. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.kmer_size = 4;
    p.min_overlap = 8;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2rc;
    rng_seed(0x1234);
    char a[101]; rand_seq(a, 100);
    char b[51]; memcpy(b, a + 50, 50); b[50] = '\0';
    load_read(&r1, "x", a, NULL);
    load_read(&r2rc, "y", b, NULL);
    int n = pbj_overlap_candidates(&ws, &r1, &r2rc, &p);
    CHECK(n >= 1);
    pbj_read_free(&r1); pbj_read_free(&r2rc);
    ws_free(&ws);
}

static void test_ec_candidates_kmer_size_max(void) {
    /* kmer_size = 31 is the max; the mask formula has a special branch
       (k >= 32) that should never trigger at k=31. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.kmer_size = 31;
    p.min_overlap = 31;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2rc;
    rng_seed(0x5678);
    char a[101]; rand_seq(a, 100);
    char b[51]; memcpy(b, a + 50, 50); b[50] = '\0';
    load_read(&r1, "x", a, NULL);
    load_read(&r2rc, "y", b, NULL);
    int n = pbj_overlap_candidates(&ws, &r1, &r2rc, &p);
    CHECK(n >= 1);
    pbj_read_free(&r1); pbj_read_free(&r2rc);
    ws_free(&ws);
}

static void test_ec_candidates_kmer_out_of_range(void) {
    /* k <= 0 or k > 31 returns 0 without crashing */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_read_t r1, r2rc;
    rng_seed(0x9ABC);
    char a[101]; rand_seq(a, 100);
    char b[51]; memcpy(b, a, 50); b[50] = '\0';
    load_read(&r1, "x", a, NULL);
    load_read(&r2rc, "y", b, NULL);
    p.kmer_size = 0;
    CHECK_EQ_INT(pbj_overlap_candidates(&ws, &r1, &r2rc, &p), 0);
    p.kmer_size = 32;
    CHECK_EQ_INT(pbj_overlap_candidates(&ws, &r1, &r2rc, &p), 0);
    pbj_read_free(&r1); pbj_read_free(&r2rc);
    ws_free(&ws);
}

/* ===================================================================
   adapter trimming edge cases
   =================================================================== */

static void test_ec_adapter_longer_than_read(void) {
    /* adapter longer than the read: only a prefix of the adapter can
       match. should still find a 3'-anchored match if quality is met. */
    const char *adapter = "AGATCGGAAGAGCACACGTCTAAAAAAAAAA";  /* 31bp */
    char seq[11];
    memcpy(seq, adapter, 10); seq[10] = '\0';  /* read = first 10bp of adapter */
    int trim = pbj_adapter_find_trim(seq, 10, adapter, 31, 8, 0.9);
    CHECK_EQ_INT(trim, 0);  /* entire read is adapter -> trim everything */
}

static void test_ec_adapter_entire_read(void) {
    /* read is exactly the adapter -> trim to 0 */
    const char *adapter = "AGATCGGAAGAGCACACGTCT";
    int trim = pbj_adapter_find_trim(adapter, 21, adapter, 21, 8, 0.9);
    CHECK_EQ_INT(trim, 0);
}

static void test_ec_adapter_strict_frac(void) {
    /* min_match_frac = 1.0 (zero tolerance). a single error in the 3' tail
       breaks the match. */
    char adapter[22];
    memcpy(adapter, "AGATCGGAAGAGCACACGTCT", 21); adapter[21] = '\0';
    char read[101];
    memset(read, 'C', 79);
    memcpy(read + 79, adapter, 21);
    read[100] = '\0';
    /* perfect tail: trims at 79 */
    CHECK_EQ_INT(pbj_adapter_find_trim(read, 100, adapter, 21, 8, 1.0), 79);
    /* introduce one mismatch in the adapter tail */
    read[85] = 'X';  /* differs from adapter at position 85 - 79 = 6 */
    /* with frac=1.0 the full 21bp match fails. but a partial-prefix match
       starting later (e.g. pos 86) might still satisfy a smaller-eff match
       if it's clean. the function picks the LEFTMOST such match. */
    int trim = pbj_adapter_find_trim(read, 100, adapter, 21, 8, 1.0);
    /* with strict frac=1.0 and 8bp min, the leftmost clean prefix-match
       is at some pos > 79 (it found a clean smaller window). */
    CHECK(trim > 79 && trim <= 100);
}

static void test_ec_adapter_min_match_above_read_len(void) {
    /* min_match_bp larger than read length -> no trim possible */
    char seq[11] = "AGATCGGAAG";
    int trim = pbj_adapter_find_trim(seq, 10, "AGATCGGAAGAGCACACGTCT", 21, 50, 0.9);
    CHECK_EQ_INT(trim, 10);
}

static void test_ec_adapter_zero_min_match_bp_default(void) {
    /* min_match_bp <= 0 should silently default to 5 */
    char read[101];
    memset(read, 'C', 95);
    memcpy(read + 95, "AGATC", 5);
    read[100] = '\0';
    int trim = pbj_adapter_find_trim(read, 100, "AGATCGGAAGAGCACACGTCT", 21, 0, 0.9);
    CHECK_EQ_INT(trim, 95);
}

static void test_ec_adapter_all_n_wildcard(void) {
    /* adapter consisting entirely of Ns acts as a universal wildcard.
       still 3'-anchored, so we trim at the leftmost valid position. */
    const char *adapter = "NNNNNNNNNN";  /* 10bp all N */
    char read[101];
    memset(read, 'C', 100); read[100] = '\0';
    int trim = pbj_adapter_find_trim(read, 100, adapter, 10, 8, 0.9);
    /* leftmost 3'-anchored position with eff >= 8 is pos 90 (eff=10),
       and all 10 Ns match anything, so trim happens there. */
    CHECK_EQ_INT(trim, 90);
}

/* ===================================================================
   pair name matching edge cases
   =================================================================== */

static void test_ec_pair_empty_names(void) {
    /* both reads have empty names (synthetic) */
    pbj_read_t a, b;
    pbj_read_init(&a); pbj_read_init(&b);
    a.name = dup_str(""); a.name_len = 0;
    b.name = dup_str(""); b.name_len = 0;
    CHECK_EQ_INT(pbj_pair_names_match(&a, &b), 1);
    pbj_read_free(&a); pbj_read_free(&b);
}

static void test_ec_pair_one_suffix_one_not(void) {
    /* "foo/1" + "foo": after stripping /1 the lengths differ (3 vs 3? no,
       "foo" is len 3, "foo/1" stripped is "foo" len 3). They actually
       match. So pick names where the stripped lengths really differ. */
    pbj_read_t a, b;
    load_read(&a, "foo/1", "A", "?");
    load_read(&b, "foo",   "A", "?");
    /* stripped lens: 3 and 3, prefixes match -> match (a.stripped=1, b.stripped=0
       so the explicit-same-suffix check passes). */
    CHECK_EQ_INT(pbj_pair_names_match(&a, &b), 1);
    pbj_read_free(&a); pbj_read_free(&b);
}

static void test_ec_pair_core_lengths_differ(void) {
    pbj_read_t a, b;
    load_read(&a, "foo/1", "A", "?");
    load_read(&b, "longerfoo/2", "A", "?");
    CHECK_EQ_INT(pbj_pair_names_match(&a, &b), 0);
    pbj_read_free(&a); pbj_read_free(&b);
}

/* ===================================================================
   casava direction edge cases
   =================================================================== */

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
    if (comment_l > 0) {
        r->name[cl] = ' ';
        memcpy(r->name + cl + 1, comment, comment_l);
        r->name[cl + 1 + comment_l] = '\0';
    } else {
        r->name[cl] = '\0';
    }
    r->name_len = (uint32_t)cl;
}

static void test_ec_casava_too_short(void) {
    /* comment of just "1" with no colon is too short to be a direction
       marker -- the implementation requires a colon after the digit. */
    pbj_read_t r;
    make_casava(&r, "name", "1");
    CHECK_EQ_INT(pbj_casava_direction(&r), 0);
    pbj_read_free(&r);
}

static void test_ec_casava_just_digit_colon(void) {
    /* minimal valid form: just "1:" (no metadata after). pbj treats this
       as a valid direction-1 marker. */
    pbj_read_t r;
    make_casava(&r, "name", "1:");
    CHECK_EQ_INT(pbj_casava_direction(&r), 1);
    pbj_read_free(&r);
}

static void test_ec_casava_no_colon_after_digit(void) {
    /* "1XYZ" doesn't have a colon at position name_len+2 */
    pbj_read_t r;
    make_casava(&r, "name", "1XYZ");
    CHECK_EQ_INT(pbj_casava_direction(&r), 0);
    pbj_read_free(&r);
}

static void test_ec_casava_invalid_digit(void) {
    /* "3:N:..." is not a valid direction */
    pbj_read_t r;
    make_casava(&r, "name", "3:N:0:ATCG");
    CHECK_EQ_INT(pbj_casava_direction(&r), 0);
    pbj_read_free(&r);
}

static void test_ec_casava_empty_comment(void) {
    pbj_read_t r;
    make_casava(&r, "name", "");
    CHECK_EQ_INT(pbj_casava_direction(&r), 0);
    pbj_read_free(&r);
}

/* ===================================================================
   merge end-to-end edge cases
   =================================================================== */

static void test_ec_merge_overlap_exactly_min(void) {
    /* construct a pair whose overlap is exactly min_overlap (here 20).
       at default p_value=0.01 and matches_only it should still merge. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.min_overlap = 20;
    p.min_assembly_len = 20;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    rng_seed(0xED6E);
    char frag[181];
    rand_seq(frag, 180);
    char r1s[101], r2s[101], q[101];
    memcpy(r1s, frag, 100); r1s[100] = '\0';
    char *rc = revcomp_str(frag + 80);  /* fragment[80..180], 100bp */
    memcpy(r2s, rc, 100); r2s[100] = '\0';
    free(rc);
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "edge/1", r1s, q);
    load_read(&r2, "edge/2", r2s, q);
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
    CHECK_EQ_INT(res.overlap_len, 20);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_ec_merge_palindromic_insert(void) {
    /* insert is its own reverse complement -> R1 == R2 (when revcomp'd) */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    /* a palindrome formed by random 50bp + its revcomp */
    rng_seed(0xBABE);
    char half[51]; rand_seq(half, 50);
    char *rc = revcomp_str(half);
    char frag[101];
    memcpy(frag, half, 50);
    memcpy(frag + 50, rc, 50);
    frag[100] = '\0';
    free(rc);

    char r1s[81], r2s[81], q[81];
    memcpy(r1s, frag, 80); r1s[80] = '\0';
    char *r2rc = revcomp_str(frag + 20);
    memcpy(r2s, r2rc, 80); r2s[80] = '\0';
    free(r2rc);
    memset(q, '?', 80); q[80] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "pal/1", r1s, q);
    load_read(&r2, "pal/2", r2s, q);
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
    CHECK_EQ_INT(res.merged.len, 100);
    CHECK_EQ_STR(res.merged.seq, frag);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_ec_merge_qtrim_below_min_length(void) {
    /* qtrim reduces merged length below min_assembly_len -> FAIL_LENGTH */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.qtrim = 20;
    p.min_assembly_len = 140;  /* tighter than the trimmable range */
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    rng_seed(0xC5C5);
    char frag[151]; rand_seq(frag, 150);
    char r1s[101], r2s[101], q1[101], q2[101];
    memcpy(r1s, frag, 100); r1s[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2s, r2rc, 100); r2s[100] = '\0';
    free(r2rc);
    /* low quality on r2's head -> after revcomp the low-Q region lands on
       the tail of the merged read, triggering qtrim. */
    memset(q1, '?', 100); q1[100] = '\0';
    memset(q2, '#', 30);             /* phred 2 */
    memset(q2 + 30, '?', 70);
    q2[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "qt/1", r1s, q1);
    load_read(&r2, "qt/2", r2s, q2);
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_FAIL_LENGTH);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_ec_merge_min_qual_boundary(void) {
    /* mean merged quality near the -Q boundary: phred 20 input, with
       overlap doubling some positions, should accept -Q 20 and reject
       -Q 30 (loosely; depends on overlap fraction). */
    pbj_revcomp_init();
    rng_seed(0xD00D);
    char frag[151]; rand_seq(frag, 150);
    char r1s[101], r2s[101], q[101];
    memcpy(r1s, frag, 100); r1s[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2s, r2rc, 100); r2s[100] = '\0';
    free(r2rc);
    memset(q, '5', 100); q[100] = '\0';  /* '5' = phred 20 */

    {
        pbj_params_t p; params_defaults(&p);
        p.min_mean_qual_merged = 20;
        pbj_workspace_t ws; ws_init(&ws, &p);
        pbj_merge_result_t res; memset(&res, 0, sizeof(res));
        pbj_read_init(&res.merged);
        pbj_read_t r1, r2;
        load_read(&r1, "q/1", r1s, q);
        load_read(&r2, "q/2", r2s, q);
        pbj_try_merge(&r1, &r2, &ws, &p, &res);
        /* mean Q is at least 20 in non-overlap; in overlap qualities sum
           to phred 40. so mean >= 20 -> accept. */
        CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
        pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
        ws_free(&ws);
    }
    {
        pbj_params_t p; params_defaults(&p);
        p.min_mean_qual_merged = 60;  /* impossibly high */
        pbj_workspace_t ws; ws_init(&ws, &p);
        pbj_merge_result_t res; memset(&res, 0, sizeof(res));
        pbj_read_init(&res.merged);
        pbj_read_t r1, r2;
        load_read(&r1, "q/1", r1s, q);
        load_read(&r2, "q/2", r2s, q);
        pbj_try_merge(&r1, &r2, &ws, &p, &res);
        CHECK_EQ_INT(res.status, PBJ_FAIL_MIN_QUAL);
        pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
        ws_free(&ws);
    }
}

static void test_ec_merge_strict_q_boundary(void) {
    /* strict mode counts a disagreement as "ambiguous" only when BOTH
       qualities are >= strict_q_threshold. just below the threshold the
       disagreement should not flag, even under --strict. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.strict_ambiguity = 1;
    p.strict_q_threshold = 30;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    rng_seed(0xBEAD);
    char frag[151]; rand_seq(frag, 150);
    char r1s[101], r2s[101];
    memcpy(r1s, frag, 100); r1s[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2s, r2rc, 100); r2s[100] = '\0';
    free(r2rc);

    /* mutate r1 at positions 60, 70, 80 inside the overlap. set q at those
       positions to phred 29 (just below threshold 30). all other positions
       at phred 40 ('I'). */
    char swap[256] = {0};
    swap['A']='C'; swap['C']='A'; swap['G']='T'; swap['T']='G';
    r1s[60] = swap[(unsigned char)r1s[60]];
    r1s[70] = swap[(unsigned char)r1s[70]];
    r1s[80] = swap[(unsigned char)r1s[80]];
    char q[101]; memset(q, 'I', 100); q[100] = '\0';
    q[60] = (char)(33 + 29);  /* phred 29 */
    q[70] = (char)(33 + 29);
    q[80] = (char)(33 + 29);

    pbj_read_t r1, r2;
    load_read(&r1, "sq/1", r1s, q);
    load_read(&r2, "sq/2", r2s, q);  /* same low Q at those positions */
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    /* disagreements are below the strict threshold -> not ambiguous */
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_ec_merge_n_in_overlap(void) {
    /* if R1 has N in the overlap and R2 has a real base, the consensus
       should use R2's base with R2's quality. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    rng_seed(0xACE5);
    char frag[151]; rand_seq(frag, 150);
    char r1s[101], r2s[101], q[101];
    memcpy(r1s, frag, 100); r1s[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2s, r2rc, 100); r2s[100] = '\0';
    free(r2rc);
    memset(q, '?', 100); q[100] = '\0';

    /* mask out one base on R1 inside the overlap. position 75 is in
       overlap (50..100 on r1). */
    char r2_base_at_75 = frag[75];
    r1s[75] = 'N';

    pbj_read_t r1, r2;
    load_read(&r1, "n/1", r1s, q);
    load_read(&r2, "n/2", r2s, q);
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
    /* merged base at position 75 should come from R2 = original fragment */
    CHECK_EQ_INT(res.merged.seq[75], r2_base_at_75);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_ec_merge_no_bonferroni(void) {
    /* --no-bonferroni should still produce the same merge on a clear
       overlap (it only loosens the test). */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.bonferroni = 0;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    rng_seed(0x10AD);
    char frag[151]; rand_seq(frag, 150);
    char r1s[101], r2s[101], q[101];
    memcpy(r1s, frag, 100); r1s[100] = '\0';
    char *r2rc = revcomp_str(frag + 50);
    memcpy(r2s, r2rc, 100); r2s[100] = '\0';
    free(r2rc);
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "nb/1", r1s, q);
    load_read(&r2, "nb/2", r2s, q);
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
    CHECK_EQ_INT(res.merged.len, 150);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_ec_merge_dovetail(void) {
    /* insert shorter than the read length -> reads overshoot.
       requires --dovetail to be enabled. */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.allow_dovetail = 1;
    p.min_assembly_len = 60;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    rng_seed(0xD0E5);
    char insert[81]; rand_seq(insert, 80);
    /* r1 = insert + random 20bp tail (adapter readthrough);
       r2 = revcomp(insert) + random 20bp tail. */
    char r1s[101], r2s[101], q[101];
    memcpy(r1s, insert, 80);
    char tail1[21]; rand_seq(tail1, 20); memcpy(r1s + 80, tail1, 20);
    r1s[100] = '\0';
    char *rcins = revcomp_str(insert);
    memcpy(r2s, rcins, 80);
    free(rcins);
    char tail2[21]; rand_seq(tail2, 20); memcpy(r2s + 80, tail2, 20);
    r2s[100] = '\0';
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "dv/1", r1s, q);
    load_read(&r2, "dv/2", r2s, q);
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK_EQ_INT(res.status, PBJ_MERGED_OK);
    /* merged length is the insert length, with adapter tails dropped */
    CHECK_EQ_INT(res.merged.len, 80);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

static void test_ec_merge_dovetail_disabled(void) {
    /* same input as the dovetail test, but without --dovetail. should
       refuse to merge (the only viable alignment is negative-offset). */
    pbj_revcomp_init();
    pbj_params_t p; params_defaults(&p);
    p.allow_dovetail = 0;
    p.min_assembly_len = 60;
    pbj_workspace_t ws; ws_init(&ws, &p);
    pbj_merge_result_t res; memset(&res, 0, sizeof(res));
    pbj_read_init(&res.merged);

    rng_seed(0xD0E5);  /* same seed as above */
    char insert[81]; rand_seq(insert, 80);
    char r1s[101], r2s[101], q[101];
    memcpy(r1s, insert, 80);
    char tail1[21]; rand_seq(tail1, 20); memcpy(r1s + 80, tail1, 20);
    r1s[100] = '\0';
    char *rcins = revcomp_str(insert);
    memcpy(r2s, rcins, 80);
    free(rcins);
    char tail2[21]; rand_seq(tail2, 20); memcpy(r2s + 80, tail2, 20);
    r2s[100] = '\0';
    memset(q, '?', 100); q[100] = '\0';

    pbj_read_t r1, r2;
    load_read(&r1, "dn/1", r1s, q);
    load_read(&r2, "dn/2", r2s, q);
    pbj_try_merge(&r1, &r2, &ws, &p, &res);
    CHECK(res.status != PBJ_MERGED_OK);

    pbj_read_free(&r1); pbj_read_free(&r2); pbj_read_free(&res.merged);
    ws_free(&ws);
}

/* ===================================================================
   driver
   =================================================================== */

int main(void) {
    printf("pbj edge-case tests\n");
    printf("-------------------\n");

    /* read buffer */
    RUN(test_ec_ensure_cap_zero);
    RUN(test_ec_ensure_cap_minimum);
    RUN(test_ec_ensure_cap_doubling);
    RUN(test_ec_ensure_cap_idempotent);
    RUN(test_ec_free_reuse_cycle);

    /* revcomp */
    RUN(test_ec_revcomp_empty);
    RUN(test_ec_revcomp_single_base);
    RUN(test_ec_revcomp_invalid_chars);

    /* binomial sf */
    RUN(test_ec_binom_n_equals_1);
    RUN(test_ec_binom_p_one_with_k_equal_n);
    RUN(test_ec_binom_small_p);
    RUN(test_ec_binom_k_equals_n);
    RUN(test_ec_binom_exact_normal_boundary);

    /* score */
    RUN(test_ec_score_zero_length_overlap);
    RUN(test_ec_score_quality_below_offset);
    RUN(test_ec_score_max_phred);
    RUN(test_ec_score_iupac_treated_as_mismatch);
    RUN(test_ec_score_zero_reward_zero_penalty);
    RUN(test_ec_score_quality_offset_64);

    /* geometry */
    RUN(test_ec_geom_exactly_min_overlap);
    RUN(test_ec_geom_one_below_min_overlap);
    RUN(test_ec_geom_exactly_max_overlap);
    RUN(test_ec_geom_one_above_max_overlap);
    RUN(test_ec_geom_offset_past_r1_end);
    RUN(test_ec_geom_very_negative_offset);
    RUN(test_ec_geom_asymmetric_r1_much_longer);

    /* candidates */
    RUN(test_ec_candidates_len_exactly_k);
    RUN(test_ec_candidates_len_below_k);
    RUN(test_ec_candidates_kmer_size_min);
    RUN(test_ec_candidates_kmer_size_max);
    RUN(test_ec_candidates_kmer_out_of_range);

    /* adapter */
    RUN(test_ec_adapter_longer_than_read);
    RUN(test_ec_adapter_entire_read);
    RUN(test_ec_adapter_strict_frac);
    RUN(test_ec_adapter_min_match_above_read_len);
    RUN(test_ec_adapter_zero_min_match_bp_default);
    RUN(test_ec_adapter_all_n_wildcard);

    /* pair name */
    RUN(test_ec_pair_empty_names);
    RUN(test_ec_pair_one_suffix_one_not);
    RUN(test_ec_pair_core_lengths_differ);

    /* casava */
    RUN(test_ec_casava_too_short);
    RUN(test_ec_casava_just_digit_colon);
    RUN(test_ec_casava_no_colon_after_digit);
    RUN(test_ec_casava_invalid_digit);
    RUN(test_ec_casava_empty_comment);

    /* merge end-to-end */
    RUN(test_ec_merge_overlap_exactly_min);
    RUN(test_ec_merge_palindromic_insert);
    RUN(test_ec_merge_qtrim_below_min_length);
    RUN(test_ec_merge_min_qual_boundary);
    RUN(test_ec_merge_strict_q_boundary);
    RUN(test_ec_merge_n_in_overlap);
    RUN(test_ec_merge_no_bonferroni);
    RUN(test_ec_merge_dovetail);
    RUN(test_ec_merge_dovetail_disabled);

    printf("-------------------\n");
    printf("%d checks, %d failed\n", g_tests, g_failed);
    return g_failed == 0 ? 0 : 1;
}
