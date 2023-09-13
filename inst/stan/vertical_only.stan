functions {
#include "include/functions.stan"
}
data {
#include "include/data_start.stan"
#include "include/data_deviation_start.stan"

  int<lower=1> N_indices;
  int<lower=1> N_blocks;
  int block_indices[N_indices];
  int block_last_index[N_blocks];
  int block_N_responses[N_blocks];

  real<lower=0> gamma_deviation_a;
  real<lower=0> gamma_deviation_b;
}
transformed data {
#include "include/transformed_data.stan"
}
parameters {
#include "include/parameters_start.stan"
#include "include/parameters_deviation_start.stan"
  vector<lower=0>[P_deviation_warping] gamma_deviation_vertical;
}
transformed parameters {
#include "include/transformed_parameters_outer_start.stan"

  {
    real x_vertical_warped[N];

#include "include/transformed_parameters_inner_start.stan"
#include "include/transformed_parameters_inner_deviation_start.stan"

    for (i in 1:N) {
      x_vertical_warped[i] = X_deviation_warping[i, :] * cumulative_sum(gamma_deviation_vertical);
    }

    for (i in 1:N_blocks) {
#include "include/transformed_parameters_inner_deviation_block_start.stan"

      K_block[:N_current_block, :N_current_block] = geowarp_process_covariance_1d(
        sigma_squared_nugget,
        deviation_sd[indices_current_block[:N_current_block]],
        x_vertical_warped[indices_current_block[:N_current_block]],
        smoothness
      );

#include "include/transformed_parameters_inner_deviation_block_end.stan"
    }

#include "include/transformed_parameters_inner_end.stan"
  }
}
model {
#include "include/model_start.stan"
#include "include/model_deviation_start.stan"

  if (gamma_deviation_b > 0) {
    gamma_deviation_vertical ~ gamma(gamma_deviation_a, gamma_deviation_b);
  }
}
