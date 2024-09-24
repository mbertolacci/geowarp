functions {
#include "include/functions_start.stan"

vector geowarp_vecchia_reduce_sum_vertical_only(
  int start,
  int end,
  int grainsize,
  real[] x_vertical_warped,
  real sigma_squared_nugget,
  vector deviation_sd,
  data vector y,
  data vector y_tilde,
  data matrix X_mean,
  data int[] block_indices,
  data int[] block_last_index,
  data int[] block_N_responses,
  data int[,] indices_X_mean_non_zero,
  data int[] P_X_mean_non_zero,
  data int N_block_max,
  data real smoothness
);

vector geowarp_vecchia_partial_sums_vertical_only(
  int start,
  int end,
  real[] x_vertical_warped,
  real sigma_squared_nugget,
  vector deviation_sd,
  data vector y,
  data vector y_tilde,
  data matrix X_mean,
  data int[] block_indices,
  data int[] block_last_index,
  data int[] block_N_responses,
  data int[,] indices_X_mean_non_zero,
  data int[] P_X_mean_non_zero,
  data int N_block_max,
  data real smoothness
) {
#include "include/functions_geowarp_vecchia_partial_sums_start.stan"

  K_block[:N_current_block, :N_current_block] = geowarp_process_covariance_1d(
    sigma_squared_nugget,
    deviation_sd[indices_current_block[:N_current_block]],
    x_vertical_warped[indices_current_block[:N_current_block]],
    smoothness
  );

#include "include/functions_geowarp_vecchia_partial_sums_end.stan"
}
}
data {
#include "include/data_start.stan"
#include "include/data_deviation_start.stan"

  int<lower=1> N_indices;
  int<lower=1> N_blocks;
  int block_indices[N_indices];
  int block_last_index[N_blocks];
  int block_N_responses[N_blocks];

  int gamma_deviation_prior_type;
  real<lower=0> gamma_deviation_lower;
  real<lower=0> gamma_deviation_upper;
  real<lower=0> gamma_deviation_a;
  real<lower=0> gamma_deviation_b;

  int use_parallel;
  int<lower=1> grain_size;
}
transformed data {
#include "include/transformed_data.stan"
}
parameters {
#include "include/parameters_start.stan"
#include "include/parameters_deviation_start.stan"
  vector<
    lower=gamma_deviation_lower,
    upper=gamma_deviation_upper
  >[P_deviation_warping] gamma_deviation_vertical;
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

    if (use_parallel > 0) {
      components = geowarp_vecchia_reduce_sum_vertical_only(
        1,
        N_blocks,
        grain_size,
        x_vertical_warped,
        sigma_squared_nugget,
        deviation_sd,
        y,
        y_tilde,
        X_mean,
        block_indices,
        block_last_index,
        block_N_responses,
        indices_X_mean_non_zero,
        P_X_mean_non_zero,
        N_block_max,
        smoothness
      );
    } else {
      components = geowarp_vecchia_partial_sums_vertical_only(
        1,
        N_blocks,
        x_vertical_warped,
        sigma_squared_nugget,
        deviation_sd,
        y,
        y_tilde,
        X_mean,
        block_indices,
        block_last_index,
        block_N_responses,
        indices_X_mean_non_zero,
        P_X_mean_non_zero,
        N_block_max,
        smoothness
      );
    }

#include "include/transformed_parameters_inner_unpack.stan"
#include "include/transformed_parameters_inner_end.stan"
  }
}
model {
#include "include/model_start.stan"
#include "include/model_deviation_start.stan"

  if (gamma_deviation_prior_type == 1) {
    gamma_deviation_vertical ~ gamma(gamma_deviation_a, gamma_deviation_b);
  } else if (gamma_deviation_prior_type == 2) {
    // Inverse uniform distribution
    target += -2.0 * sum(log(gamma_deviation_vertical));
  } else {
    gamma_deviation_vertical ~ uniform(gamma_deviation_lower, gamma_deviation_upper);
  }
}
