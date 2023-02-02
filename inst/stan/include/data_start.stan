int<lower=1> N;
int<lower=2> D;
vector[D] x[N];
vector[N] y;
vector[D] scaling;

int<lower=1> P_mean_fixed;
matrix[N, P_mean_fixed] X_mean_fixed;
vector[P_mean_fixed] alpha_mean;
matrix[P_mean_fixed, P_mean_fixed] alpha_precision;

int<lower=0> P_mean_random;
matrix[N, P_mean_random] X_mean_random;
real<lower=0> tau_squared_mean_random_a;
real<lower=0> tau_squared_mean_random_b;

real<lower=0> sigma_squared_nugget_a;
real<lower=0> sigma_squared_nugget_b;
