int<lower=1> P_deviation_warping;
matrix[N, P_deviation_warping] X_deviation_warping;

int<lower=1> P_deviation_fixed;
matrix[N, P_deviation_fixed] X_deviation_fixed;
vector[P_deviation_fixed] eta_deviation_mean;
cov_matrix[P_deviation_fixed] eta_deviation_precision;

int<lower=0> P_deviation_random;
matrix[N, P_deviation_random] X_deviation_random;
real<lower=0> delta_deviation_random;
real<lower=0> ell_deviation_random_scale;
real<lower=0> tau_squared_deviation_random_scale;
