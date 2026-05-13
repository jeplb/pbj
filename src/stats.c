#include "stats.h"
#include <math.h>

/* exact binomial survival function: p(x >= k) for x ~ binomial(n, p).
   computed in log-space using lgamma to stay stable for n up to ~256.
   for n above that, normal approximation is fine and the test is
   effectively certain anyway. */
static double binomial_sf_exact(int k, int n, double p) {
    if (k <= 0) return 1.0;
    if (k > n) return 0.0;
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1.0;

    double log_p   = log(p);
    double log_1mp = log(1.0 - p);
    double log_n_fac = lgamma((double)n + 1.0);

    /* sum from k to n. start from k and accumulate. */
    double sum = 0.0;
    for (int i = k; i <= n; i++) {
        double log_choose = log_n_fac
                          - lgamma((double)i + 1.0)
                          - lgamma((double)(n - i) + 1.0);
        double log_term = log_choose + (double)i * log_p + (double)(n - i) * log_1mp;
        sum += exp(log_term);
    }
    if (sum > 1.0) sum = 1.0;
    return sum;
}

static double normal_sf(double z) {
    /* p(z >= z) = 0.5 * erfc(z / sqrt(2)) */
    return 0.5 * erfc(z / sqrt(2.0));
}

double pbj_binomial_sf(int k, int n, double p) {
    if (n <= 0) return 1.0;
    if (n <= 256) return binomial_sf_exact(k, n, p);

    double mu = (double)n * p;
    double sd = sqrt((double)n * p * (1.0 - p));
    if (sd <= 0.0) return (k > mu) ? 0.0 : 1.0;
    /* continuity correction */
    double z = ((double)k - 0.5 - mu) / sd;
    return normal_sf(z);
}
