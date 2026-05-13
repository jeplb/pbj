#ifndef PBJ_STATS_H
#define PBJ_STATS_H

/* upper-tail binomial p-value: p(x >= k | n, p_random_match).
   uses exact summation for n <= 256, normal approximation otherwise. */
double pbj_binomial_sf(int k, int n, double p);

#endif
