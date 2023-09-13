// Body in stan_meta_header.hpp
matrix geowarp_process_covariance(
  real sigma_squared_nugget,
  vector deviation_sd,
  vector[] x,
  data real smoothness
);

// Body in stan_meta_header.hpp
matrix geowarp_process_covariance_1d(
  real sigma_squared_nugget,
  vector deviation_sd,
  real[] x,
  data real smoothness
);

// Body in stan_meta_header.hpp
int get_N_block_max(int[] block_last_index);

matrix rw1d_precision(int n, real sigma_squared) {
  matrix[n, n] output = diag_matrix(rep_vector(2, n));
  if (n == 0) {
    return output;
  }
  for (i in 2:n) {
    output[i, i - 1] = -1;
    output[i - 1, i] = -1;
  }
  output[n, n] = 1;
  return output / sigma_squared;
}

matrix exp1d_precision(int n, real delta, real ell, real sigma_squared) {
  matrix[n, n] output = rep_matrix(0, n, n);
  real e_lambda_d_x;
  real e_2lambda_d_x;
  real r_major;
  real r_minor;

  if (n == 1) {
    output[1, 1] = 1 / sigma_squared;
  }
  if (n <= 1) {
    return output;
  }
  e_lambda_d_x = exp(-delta / ell);
  e_2lambda_d_x = exp(-2 * delta / ell);
  r_major = e_2lambda_d_x / (1 - e_2lambda_d_x);
  r_minor = e_lambda_d_x / (1 - e_2lambda_d_x);

  output[1, 1] = 1 + r_major;
  for (i in 2:(n-1)) {
    output[i, i] = 1 + 2 * r_major;
    output[i, i - 1] = -r_minor;
    output[i - 1, i] = -r_minor;
  }
  output[n, n - 1] = -r_minor;
  output[n - 1, n] = -r_minor;
  output[n, n] = 1 + r_major;
  return output / sigma_squared;
}

// Solves (LL')X = b for X, via L^{-T} L^{-1} b
vector chol_solve_L_b(
  matrix L,
  vector b
) {
  return mdivide_right_tri_low(mdivide_left_tri_low(L, b)', L)';
}
