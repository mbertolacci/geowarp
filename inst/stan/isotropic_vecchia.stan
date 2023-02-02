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

  vector<lower=0>[D] gamma_deviation_a;
  vector<lower=0>[D] gamma_deviation_b;
}
transformed data {
#include "include/transformed_data_vecchia.stan"
}
parameters {
#include "include/parameters_start.stan"
#include "include/parameters_deviation_start.stan"
  vector<lower=0>[D - 1] gamma_deviation_horizontal;
  vector<lower=0>[P_deviation_warping] gamma_deviation_vertical;
}
transformed parameters {
#include "include/transformed_parameters_outer_start.stan"

  {
    vector[D] x_warped[N];

#include "include/transformed_parameters_inner_start.stan"
#include "include/transformed_parameters_inner_deviation_vecchia_start.stan"

    for (i in 1:N) {
      for (j in 1:(D-1)) {
        x_warped[i][j] = gamma_deviation_horizontal[j] * scaling[j] * x[i][j];
      }
      x_warped[i][D] = X_deviation_warping[i, :] * cumulative_sum(gamma_deviation_vertical);
    }

    for (i in 1:N_blocks) {
#include "include/transformed_parameters_inner_deviation_vecchia_block_start.stan"

      K_block[:N_current_block, :N_current_block] = exponential_cov_heteroskedastic_vector(
        x_warped[indices_current_block[:N_current_block]],
        deviation_sd[indices_current_block[:N_current_block]],
        sigma_squared_nugget
      );

#include "include/transformed_parameters_inner_deviation_vecchia_block_end.stan"
    }

#include "include/transformed_parameters_inner_end.stan"
  }
}
model {
#include "include/model_start.stan"
#include "include/model_deviation_start.stan"

  for (i in 1:(D-1)) {
    if (gamma_deviation_b[i] > 0) {
      gamma_deviation_horizontal[i] ~ gamma(gamma_deviation_a[i], gamma_deviation_b[i]);
    }
  }
  if (gamma_deviation_b[D] > 0) {
    gamma_deviation_vertical ~ gamma(gamma_deviation_a[D], gamma_deviation_b[D]);
  }
}
generated quantities {
#include "include/generated_quantities.stan"
}
