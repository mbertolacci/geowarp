// matrix geowarp_process_covariance_naive(
//   real sigma_squared_nugget,
//   vector deviation_sd,
//   vector[] x,
//   data real smoothness
// ) {
//   int n = size(x);
//   matrix[n, n] output;
//   if (smoothness == 0.5) {
//     for (j in 1:n) {
//       output[j, j] = square(deviation_sd[j]) + sigma_squared_nugget;
//       for (k in (j+1):n) {
//         output[j, k] = (
//           deviation_sd[j]
//            * deviation_sd[k]
//            * exp(-distance(x[j], x[k]))
//         );
//         output[k, j] = output[j, k];
//       }
//     }
//   } else if (smoothness == 1.5) {
//     for (j in 1:n) {
//       output[j, j] = square(deviation_sd[j]) + sigma_squared_nugget;
//       for (k in (j+1):n) {
//         real d = distance(x[j], x[k]);
//         output[j, k] = (
//           deviation_sd[j]
//            * deviation_sd[k]
//            * (1.0 + sqrt(3) * d)
//            * exp(-sqrt(3) * d)
//         );
//         output[k, j] = output[j, k];
//       }
//     }
//   }
//   return output;
// }

// matrix geowarp_process_covariance_1d_naive(
//   real sigma_squared_nugget,
//   vector deviation_sd,
//   real[] x,
//   data real smoothness
// ) {
//   int n = size(x);
//   matrix[n, n] output;
//   if (smoothness == 0.5) {
//     for (j in 1:n) {
//       output[j, j] = square(deviation_sd[j]) + sigma_squared_nugget;
//       for (k in (j+1):n) {
//         output[j, k] = (
//           deviation_sd[j]
//            * deviation_sd[k]
//            * exp(-fabs(x[j] - x[k]))
//         );
//         output[k, j] = output[j, k];
//       }
//     }
//   } else if (smoothness == 1.5) {
//     for (j in 1:n) {
//       output[j, j] = square(deviation_sd[j]) + sigma_squared_nugget;
//       for (k in (j+1):n) {
//         real d = fabs(x[j] - x[k]);
//         output[j, k] = (
//           deviation_sd[j]
//            * deviation_sd[k]
//            * (1.0 + sqrt(3) * d)
//            * exp(-sqrt(3) * d)
//         );
//         output[k, j] = output[j, k];
//       }
//     }
//   }
//   return output;
// }

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

int get_block_start(int i, int[] block_last_index) {
  if (i == 1) {
    return 1;
  } else {
    return block_last_index[i - 1] + 1;
  }
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
