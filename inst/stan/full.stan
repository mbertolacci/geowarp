functions {
#include "include/functions.stan"
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

    for (i in 1:N_blocks) {
#include "include/transformed_parameters_inner_deviation_block_start.stan"

      K_block[:N_current_block, :N_current_block] = geowarp_process_covariance(
        sigma_squared_nugget,
        deviation_sd[indices_current_block[:N_current_block]],
        x_warped[indices_current_block[:N_current_block]],
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
