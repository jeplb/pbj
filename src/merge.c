#include "merge.h"
#include "overlap.h"
#include "score.h"
#include "stats.h"
#include "revcomp.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define MAX_PHRED 93

/* posterior probability that the picked (higher-q) base is the true base
   given that the two reads disagree, under independent-error model and
   uniform prior over the four bases. derived from bayes rule with
   e1 = 10^(-q_picked/10), e2 = 10^(-q_other/10):

     p(true = picked) is proportional to (1 - e1) * (e2 / 3)
     p(true = other ) is proportional to (e1 / 3) * (1 - e2)
     p(true = third ) is proportional to (e1 / 3) * (e2 / 3), times 2 for the two such bases

   normalize and convert (1 - p_picked) back to phred. closed-form, no
   special functions beyond pow / log10. */
static int bayesian_mismatch_quality(int q_picked, int q_other) {
    if (q_picked <= 0 && q_other <= 0) return 0;
    double e1 = pow(10.0, -(double)q_picked / 10.0);
    double e2 = pow(10.0, -(double)q_other  / 10.0);
    double w_picked = (1.0 - e1) * (e2 / 3.0);
    double w_other  = (e1 / 3.0) * (1.0 - e2);
    double w_third  = (e1 / 3.0) * (e2 / 3.0);
    double denom = w_picked + w_other + 2.0 * w_third;
    if (denom <= 0.0) return 0;
    double p_correct = w_picked / denom;
    double err = 1.0 - p_correct;
    if (err <= 0.0) return MAX_PHRED;
    if (err >= 1.0) return 0;
    int q = (int)(-10.0 * log10(err) + 0.5);
    if (q < 0)         q = 0;
    if (q > MAX_PHRED) q = MAX_PHRED;
    return q;
}

/* copy the core id (kseq name field, no comment) from src into dst,
   stripping a trailing /1 or /2 suffix if present. */
static void copy_core_name(pbj_read_t *dst, const pbj_read_t *src) {
    if (!src->name || src->name_len == 0) {
        free(dst->name);
        dst->name = NULL;
        dst->name_len = 0;
        return;
    }
    uint32_t L = src->name_len;
    if (L >= 2 && src->name[L - 2] == '/' &&
        (src->name[L - 1] == '1' || src->name[L - 1] == '2')) {
        L -= 2;
    }
    char *new_name = (char*)realloc(dst->name, L + 1);
    if (!new_name) {
        fprintf(stderr, "pbj: out of memory in copy_core_name\n");
        abort();
    }
    dst->name = new_name;
    memcpy(dst->name, src->name, L);
    dst->name[L] = '\0';
    dst->name_len = L;
}

void pbj_try_merge(const pbj_read_t *r1,
                   const pbj_read_t *r2,
                   pbj_workspace_t *ws,
                   const pbj_params_t *p,
                   pbj_merge_result_t *result) {
    result->status          = PBJ_NO_CANDIDATE;
    result->matches         = 0;
    result->mismatches      = 0;
    result->overlap_len     = 0;
    result->offset          = 0;
    result->n_candidates    = 0;
    result->score           = 0.0;
    result->p_value         = 1.0;
    result->expected_errors = 0.0;

    pbj_revcomp(r2, &ws->r2_rc);

    int n_cand = pbj_overlap_candidates(ws, r1, &ws->r2_rc, p);
    result->n_candidates = n_cand;
    if (n_cand == 0) {
        result->status = PBJ_NO_CANDIDATE;
        return;
    }

    pbj_score_t best_sc = {0};
    int best_s = 0, best_r1s = 0, best_r2s = 0, best_olen = 0;
    double best_value = -1e30;
    int found = 0;

    for (int i = 0; i < n_cand; i++) {
        int s = ws->candidates[i];
        int r1s, r2s, olen;
        if (!pbj_overlap_geometry(s, (int)r1->len, (int)ws->r2_rc.len,
                                   p->min_overlap, p->max_overlap,
                                   &r1s, &r2s, &olen)) continue;
        pbj_score_t sc;
        pbj_score_overlap(r1, r1s, &ws->r2_rc, r2s, olen, p, &sc);
        if (sc.score > best_value) {
            best_value = sc.score;
            best_sc    = sc;
            best_s     = s;
            best_r1s   = r1s;
            best_r2s   = r2s;
            best_olen  = olen;
            found      = 1;
        }
    }

    if (!found) {
        result->status = PBJ_NO_CANDIDATE;
        return;
    }

    int merged_len;
    if (best_s >= 0) {
        merged_len = best_s + (int)ws->r2_rc.len;
    } else {
        merged_len = best_olen;
    }
    if (merged_len < p->min_assembly_len ||
        (p->max_assembly_len > 0 && merged_len > p->max_assembly_len)) {
        result->status = PBJ_FAIL_LENGTH;
        result->offset = best_s;
        result->overlap_len = best_olen;
        result->matches = best_sc.matches;
        result->mismatches = best_sc.mismatches;
        result->score = best_value;
        return;
    }

    /* seed-adjusted binomial test. condition on the kmer seed providing k
       guaranteed matches; the remaining (useful - k) free positions are
       i.i.d. bernoulli(p_match) under null. apply a bonferroni correction
       across the n_cand candidates we screened. */
    int useful = best_sc.matches + best_sc.mismatches;
    int n_free = useful - p->kmer_size;
    int m_free = best_sc.matches - p->kmer_size;
    if (n_free < 0) n_free = 0;
    if (m_free < 0) m_free = 0;
    double p_match = p->p_match > 0.0 ? p->p_match : 0.25;
    double p_val = (n_free > 0) ? pbj_binomial_sf(m_free, n_free, p_match) : 1.0;
    if (p->bonferroni && n_cand > 1) {
        p_val *= (double)n_cand;
        if (p_val > 1.0) p_val = 1.0;
    }
    if (p_val > p->p_value) {
        result->status = PBJ_FAIL_PVALUE;
        result->p_value = p_val;
        result->offset = best_s;
        result->overlap_len = best_olen;
        result->matches = best_sc.matches;
        result->mismatches = best_sc.mismatches;
        result->score = best_value;
        return;
    }

    pbj_read_ensure_cap(&result->merged, (uint32_t)merged_len + 1);
    result->merged.len = (uint32_t)merged_len;

    char *out_s = result->merged.seq;
    char *out_q = result->merged.qual;
    int idx     = 0;
    int ambig   = 0;
    int64_t sum_q = 0;
    int off     = p->quality_offset;
    int strict_q = p->strict_q_threshold > 0 ? p->strict_q_threshold : 20;

    if (best_s > 0) {
        memcpy(out_s + idx, r1->seq, (size_t)best_s);
        /* output qual always in phred+33; translate from input offset */
        for (int i = 0; i < best_s; i++) {
            int qv = (int)(unsigned char)r1->qual[i] - off;
            if (qv < 0)         qv = 0;
            if (qv > MAX_PHRED) qv = MAX_PHRED;
            out_q[idx + i] = (char)(qv + 33);
            sum_q += qv;
        }
        idx += best_s;
    }

    {
        const char *s1 = r1->seq  + best_r1s;
        const char *q1 = r1->qual + best_r1s;
        const char *s2 = ws->r2_rc.seq  + best_r2s;
        const char *q2 = ws->r2_rc.qual + best_r2s;
        for (int i = 0; i < best_olen; i++) {
            char b1 = s1[i];
            char b2 = s2[i];
            int qv1 = (int)(unsigned char)q1[i] - off; if (qv1 < 0) qv1 = 0;
            int qv2 = (int)(unsigned char)q2[i] - off; if (qv2 < 0) qv2 = 0;

            char out_b;
            int  out_qv;
            int  n1 = (b1 == 'N');
            int  n2 = (b2 == 'N');

            if (n1 && n2)      { out_b = 'N'; out_qv = 0; }
            else if (n1)       { out_b = b2;  out_qv = qv2; }
            else if (n2)       { out_b = b1;  out_qv = qv1; }
            else if (b1 == b2) {
                out_b  = b1;
                out_qv = qv1 + qv2;
                if (out_qv > MAX_PHRED) out_qv = MAX_PHRED;
            } else {
                int qp, qo;
                if (qv1 >= qv2) { out_b = b1; qp = qv1; qo = qv2; }
                else            { out_b = b2; qp = qv2; qo = qv1; }
                out_qv = bayesian_mismatch_quality(qp, qo);
                if (qv1 >= strict_q && qv2 >= strict_q) ambig++;
            }
            out_s[idx] = out_b;
            out_q[idx] = (char)(out_qv + 33);
            sum_q += out_qv;
            idx++;
        }
    }

    if (p->strict_ambiguity && ambig > 0) {
        result->status = PBJ_FAIL_AMBIGUOUS;
        result->offset = best_s;
        result->overlap_len = best_olen;
        result->matches = best_sc.matches;
        result->mismatches = best_sc.mismatches;
        result->score = best_value;
        result->p_value = p_val;
        return;
    }

    if (best_s >= 0) {
        int suf_start = best_r2s + best_olen;
        int suf_len   = (int)ws->r2_rc.len - suf_start;
        if (suf_len > 0) {
            memcpy(out_s + idx, ws->r2_rc.seq + suf_start, (size_t)suf_len);
            for (int i = 0; i < suf_len; i++) {
                int qv = (int)(unsigned char)ws->r2_rc.qual[suf_start + i] - off;
                if (qv < 0)         qv = 0;
                if (qv > MAX_PHRED) qv = MAX_PHRED;
                out_q[idx + i] = (char)(qv + 33);
                sum_q += qv;
            }
            idx += suf_len;
        }
    }

    out_s[idx] = '\0';
    out_q[idx] = '\0';

    /* post-merge 3' quality trim. trim while the trailing base's phred is
       below qtrim threshold. operates on +33-encoded output qualities.
       a zero-length result is always treated as fail_length so we never
       emit an empty fastq record, even when min_assembly_len is 0. */
    if (p->qtrim > 0) {
        int qt = p->qtrim;
        while (merged_len > 0 && ((int)(unsigned char)out_q[merged_len - 1] - 33) < qt) {
            sum_q -= (int)(unsigned char)out_q[merged_len - 1] - 33;
            merged_len--;
        }
        result->merged.len = (uint32_t)merged_len;
        out_s[merged_len] = '\0';
        out_q[merged_len] = '\0';
        if (merged_len == 0 || merged_len < p->min_assembly_len) {
            result->status = PBJ_FAIL_LENGTH;
            result->offset = best_s;
            result->overlap_len = best_olen;
            result->matches = best_sc.matches;
            result->mismatches = best_sc.mismatches;
            result->score = best_value;
            result->p_value = p_val;
            return;
        }
    }

    if (p->min_mean_qual_merged > 0 && merged_len > 0) {
        double mean_q = (double)sum_q / (double)merged_len;
        if (mean_q < (double)p->min_mean_qual_merged) {
            result->status = PBJ_FAIL_MIN_QUAL;
            result->offset = best_s;
            result->overlap_len = best_olen;
            result->matches = best_sc.matches;
            result->mismatches = best_sc.mismatches;
            result->score = best_value;
            result->p_value = p_val;
            return;
        }
    }

    /* expected-error filter (usearch/vsearch-style):
       ee = sum over positions of 10^(-phred/10). */
    if (p->max_ee > 0.0 && merged_len > 0) {
        double ee = 0.0;
        for (int i = 0; i < merged_len; i++) {
            int q = (int)(unsigned char)out_q[i] - 33;
            if (q < 0) q = 0;
            ee += pow(10.0, -(double)q / 10.0);
        }
        result->expected_errors = ee;
        if (ee > p->max_ee) {
            result->status = PBJ_FAIL_MAX_EE;
            result->offset = best_s;
            result->overlap_len = best_olen;
            result->matches = best_sc.matches;
            result->mismatches = best_sc.mismatches;
            result->score = best_value;
            result->p_value = p_val;
            return;
        }
    }

    /* low-complexity filter: reject if any single base (including n)
       exceeds the fraction threshold. catches homopolymers, n-rich
       reads, and very repetitive sequences. */
    if (p->max_mono_frac > 0.0 && merged_len > 0) {
        int counts[5] = {0, 0, 0, 0, 0};   /* a, c, g, t, n */
        for (int i = 0; i < merged_len; i++) {
            char b = out_s[i];
            switch (b) {
                case 'A': counts[0]++; break;
                case 'C': counts[1]++; break;
                case 'G': counts[2]++; break;
                case 'T': counts[3]++; break;
                case 'N': counts[4]++; break;
                default:  counts[4]++; break;   /* iupac codes count toward n */
            }
        }
        int maxc = counts[0];
        for (int i = 1; i < 5; i++) if (counts[i] > maxc) maxc = counts[i];
        if ((double)maxc / (double)merged_len > p->max_mono_frac) {
            result->status = PBJ_FAIL_MONO;
            result->offset = best_s;
            result->overlap_len = best_olen;
            result->matches = best_sc.matches;
            result->mismatches = best_sc.mismatches;
            result->score = best_value;
            result->p_value = p_val;
            return;
        }
    }

    copy_core_name(&result->merged, r1);

    result->status      = PBJ_MERGED_OK;
    result->offset      = best_s;
    result->overlap_len = best_olen;
    result->matches     = best_sc.matches;
    result->mismatches  = best_sc.mismatches;
    result->score       = best_value;
    result->p_value     = p_val;
}
