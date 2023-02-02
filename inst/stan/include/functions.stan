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

matrix rw1d_precision_varying(vector sigma_squared) {
  int n = rows(sigma_squared);
  matrix[n, n] output = rep_matrix(0, n, n);
  for (i in 1:(n-1)) {
    output[i, i] = (sigma_squared[i] + sigma_squared[i + 1]) / (sigma_squared[i] * sigma_squared[i + 1]);
    output[i + 1, i] = -1 / sigma_squared[i + 1];
    output[i, i + 1] = -1 / sigma_squared[i + 1];
  }
  output[n, n] = 1 / sigma_squared[n];
  return output;
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

matrix exp1d_covariance(int n, real delta, real ell, real tau_squared, real sigma_squared) {
  matrix[n, n] output = rep_matrix(0, n, n);
  for (i in 1:n) {
    for (j in (i+1):n) {
      output[i, j] = tau_squared * exp(-(
        delta * (j - i)
      ) / ell);
      output[j, i] = output[i, j];
    }
    output[i, i] += tau_squared + sigma_squared;
  }
  return output;
}

// Solves (LL')X = b for X, via L^{-T} L^{-1} b
vector chol_solve_L_b(
  matrix L,
  vector b
) {
  return mdivide_right_tri_low(mdivide_left_tri_low(L, b)', L)';
}

row_vector n_std_normals_rng(int n) {
  row_vector[n] z;
  for (i in 1:n) {
    z[i] = normal_rng(0, 1);
  }
  return z;
}

matrix exponential_cov_heteroskedastic_scalar(
  real[] x,
  vector standard_deviations,
  real sigma_squared_nugget
) {
  int n = size(x);
  matrix[n, n] output;
  for (j in 1:n) {
    output[j, j] = square(standard_deviations[j]) + sigma_squared_nugget;
    for (k in (j+1):n) {
      output[j, k] = (
        standard_deviations[j]
         * standard_deviations[k]
         * exp(-fabs(x[j] - x[k]))
      );
      output[k, j] = output[j, k];
    }
  }
  return output;
}

matrix exponential_cov_heteroskedastic_vector(
  vector[] x,
  vector standard_deviations,
  real sigma_squared_nugget
) {
  int n = size(x);
  matrix[n, n] output;
  for (j in 1:n) {
    output[j, j] = square(standard_deviations[j]) + sigma_squared_nugget;
    for (k in (j+1):n) {
      output[j, k] = (
        standard_deviations[j]
         * standard_deviations[k]
         * exp(-distance(x[j], x[k]))
      );
      output[k, j] = output[j, k];
    }
  }
  return output;
}

int get_N_block_max(int[] block_last_index) {
  int N_blocks = size(block_last_index);
  int current_block_start_d = 1;
  int N_block_max = 0;
  for (i in 1:N_blocks) {
    int N_block_i = block_last_index[i] - current_block_start_d + 1;
    if (N_block_i > N_block_max) {
      N_block_max = N_block_i;
    }
    current_block_start_d = block_last_index[i] + 1;
  }
  return N_block_max;
}
