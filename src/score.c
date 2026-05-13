#include "score.h"

void pbj_score_overlap(const pbj_read_t *r1, int r1_start,
                        const pbj_read_t *r2_rc, int r2_start,
                        int overlap_len,
                        const pbj_params_t *p,
                        pbj_score_t *out) {
    int matches    = 0;
    int mismatches = 0;
    int n_count    = 0;
    double score   = 0.0;

    int off = p->quality_offset;
    int mr  = p->match_reward;
    int mp  = p->mismatch_penalty;

    const char *s1 = r1->seq + r1_start;
    const char *q1 = r1->qual + r1_start;
    const char *s2 = r2_rc->seq + r2_start;
    const char *q2 = r2_rc->qual + r2_start;

    for (int i = 0; i < overlap_len; i++) {
        char b1 = s1[i];
        char b2 = s2[i];
        if (b1 == 'N' || b1 == 'n' || b2 == 'N' || b2 == 'n') {
            n_count++;
            continue;
        }
        int qa = (int)(unsigned char)q1[i] - off; if (qa < 0) qa = 0;
        int qb = (int)(unsigned char)q2[i] - off; if (qb < 0) qb = 0;
        int q  = (qa < qb) ? qa : qb;
        if (b1 == b2) {
            matches++;
            score += (double)mr * q;
        } else {
            mismatches++;
            score -= (double)mp * q;
        }
    }
    out->matches    = matches;
    out->mismatches = mismatches;
    out->n_count    = n_count;
    out->score      = score;
}
