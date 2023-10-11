functions {
#include "include/functions_start.stan"

vector geowarp_vecchia_reduce_sum_full(
  int start,
  int end,
  int grainsize,
  vector[] x_warped,
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

vector geowarp_vecchia_partial_sums_full(
  int start,
  int end,
  vector[] x_warped,
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

    K_block = geowarp_process_covariance(
      sigma_squared_nugget,
      deviation_sd[indices_current_block],
      x_warped[indices_current_block],
      smoothness
    );

  #include "include/functions_geowarp_vecchia_partial_sums_end.stan"
}
}
data {
#include "include/data_start.stan"
#include "include/data_deviation_start.stan"
  int<lower=1> D_horizontal_warpings;
  int<lower=0> D_geometric;
  int<lower=1,upper=D_horizontal_warpings + 1> axial_warping_unit_mapping[D];

  int<lower=1> N_indices;
  int<lower=1> N_blocks;
  int block_indices[N_indices];
  int block_last_index[N_blocks];
  int block_N_responses[N_blocks];

  vector<lower=0>[D_horizontal_warpings + 1] gamma_deviation_a;
  vector<lower=0>[D_horizontal_warpings + 1] gamma_deviation_b;
  real<lower=0> L_deviation_shape;

  int use_parallel;
  int<lower=1> grain_size;
}
transformed data {
#include "include/transformed_data.stan"
}
parameters {
#include "include/parameters_start.stan"
#include "include/parameters_deviation_start.stan"
  vector<lower=0>[D_horizontal_warpings] gamma_deviation_horizontal;
  vector<lower=0>[P_deviation_warping] gamma_deviation_vertical;
  cholesky_factor_corr[D_geometric] L_deviation;
}
transformed parameters {
#include "include/transformed_parameters_outer_start.stan"

  {
    vector[D] x_warped[N];

#include "include/transformed_parameters_inner_start.stan"
#include "include/transformed_parameters_inner_deviation_start.stan"

    for (i in 1:N) {
      for (j in 1:(D-1)) {
        x_warped[i][j] = gamma_deviation_horizontal[
          axial_warping_unit_mapping[j]
        ] * scaling[j] * x[i][j];
      }
      x_warped[i][D] = X_deviation_warping[i, :] * cumulative_sum(gamma_deviation_vertical);
      if (D_geometric > 0) {
        x_warped[i] = L_deviation' * x_warped[i];
      }
    }

    if (use_parallel > 0) {
      components = geowarp_vecchia_reduce_sum_full(
        1,
        N_blocks,
        grain_size,
        x_warped,
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
      components = geowarp_vecchia_partial_sums_full(
        1,
        N_blocks,
        x_warped,
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

  for (i in 1:D_horizontal_warpings) {
    if (gamma_deviation_b[i] > 0) {
      gamma_deviation_horizontal[i] ~ gamma(gamma_deviation_a[i], gamma_deviation_b[i]);
    }
  }
  if (gamma_deviation_b[D_horizontal_warpings + 1] > 0) {
    gamma_deviation_vertical ~ gamma(gamma_deviation_a[D_horizontal_warpings + 1], gamma_deviation_b[D_horizontal_warpings + 1]);
  }
  if (D_geometric > 0) {
    L_deviation ~ lkj_corr_cholesky(L_deviation_shape);
  }
}
