functions {
#include "include/functions.stan"
 }
data {
  int<lower=1> N;
  int<lower=2> D;
  vector[D] x[N];
  vector[N] y;

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

  int<lower=1> P_deviation_fixed;
  matrix[N, P_deviation_fixed] X_deviation_fixed;
  vector[P_deviation_fixed] eta_deviation_mean;
  cov_matrix[P_deviation_fixed] eta_deviation_precision;

  int<lower=0> P_deviation_random;
  matrix[N, P_deviation_random] X_deviation_random;
  real<lower=0> delta_deviation_random;
  real<lower=0> ell_deviation_random_scale;
  real<lower=0> tau_squared_deviation_random_a;
  real<lower=0> tau_squared_deviation_random_b;

  int<lower=1> N_indices;
  int<lower=1> N_blocks;
  int block_indices[N_indices];
  int block_last_index[N_blocks];
  int block_N_responses[N_blocks];
}
transformed data {
#include "include/transformed_data.stan"
}
parameters {
#include "include/parameters_start.stan"
#include "include/parameters_deviation_start.stan"
}
transformed parameters {
#include "include/transformed_parameters_outer_start.stan"

  {
#include "include/transformed_parameters_inner_start.stan"
#include "include/transformed_parameters_inner_deviation_start.stan"

    for (i in 1:N_blocks) {
      N_current_block = block_last_index[i] - current_block_start + 1;
      N_parents_current_block = N_current_block - block_N_responses[i];

      indices_current_block[:N_current_block] = block_indices[
        current_block_start:block_last_index[i]
      ];

      for (j in 1:P_X_mean_non_zero[i]) {
        L_inv_X_block[:N_current_block, j] = X_mean[
          indices_current_block[:N_current_block],
          indices_X_mean_non_zero[j, i]
        ] ./ deviation_sd[indices_current_block[:N_current_block]];
      }
      L_inv_y_block[:N_current_block] = (
        y[indices_current_block[:N_current_block]]
        ./ deviation_sd[indices_current_block[:N_current_block]]
      );
      L_inv_y_tilde_block[:N_current_block] = (
        y_tilde[indices_current_block[:N_current_block]]
        ./ deviation_sd[indices_current_block[:N_current_block]]
      );

      log_det += (
        2 * sum(log(deviation_sd[indices_current_block[:N_current_block]]))
      );
      y_tildet_Q_y_tilde += (
        sum(square(L_inv_y_tilde_block[(N_parents_current_block+1):N_current_block]))
      );
      Xt_Q_y[indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]] += (
        L_inv_X_block[(N_parents_current_block+1):N_current_block, :P_X_mean_non_zero[i]]'
        * L_inv_y_block[(N_parents_current_block+1):N_current_block]
      );
      Xt_Q_y_tilde[indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]] += (
        L_inv_X_block[(N_parents_current_block+1):N_current_block, :P_X_mean_non_zero[i]]'
        * L_inv_y_tilde_block[(N_parents_current_block+1):N_current_block]
      );
      Xt_Q_X[
        indices_X_mean_non_zero[:P_X_mean_non_zero[i], i],
        indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]
      ] += (
        crossprod(L_inv_X_block[
          (N_parents_current_block+1):N_current_block,
          :P_X_mean_non_zero[i]
        ])
      );

      current_block_start = block_last_index[i] + 1;
    }

#include "include/transformed_parameters_inner_end.stan"
  }
}
model {
#include "include/model_start.stan"
#include "include/model_deviation_start.stan"
}
