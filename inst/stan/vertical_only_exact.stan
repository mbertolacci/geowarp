functions {
#include "include/functions.stan"
}
data {
#include "include/data_start.stan"
#include "include/data_deviation_start.stan"

  int<lower=1> N_individuals;
  int start[N_individuals];

  real<lower=0> gamma_deviation_a;
  real<lower=0> gamma_deviation_b;
}
transformed data {
#include "include/transformed_data_exact.stan"
}
parameters {
#include "include/parameters_start.stan"
#include "include/parameters_deviation_start.stan"
  vector<lower=0>[P_deviation_warping] gamma_deviation_vertical;
}
transformed parameters {
#include "include/transformed_parameters_outer_start.stan"
  {
    matrix[N, N] K = rep_matrix(0, N, N);
    real x_vertical_warped[N];

#include "include/transformed_parameters_inner_start.stan"
#include "include/transformed_parameters_inner_deviation_exact_start.stan"

    for (i in 1:N) {
      x_vertical_warped[i] = X_deviation_warping[i, :] * cumulative_sum(gamma_deviation_vertical);
    }

    for (i in 1:N_individuals) {
      int end;
      if (i == N_individuals) {
        end = N;
      } else {
        end = start[i + 1] - 1;
      }
      K[
        start[i]:end,
        start[i]:end
      ] = exponential_cov_heteroskedastic_scalar(
        x_vertical_warped[start[i]:end],
        deviation_sd[start[i]:end],
        sigma_squared_nugget
      );
    }

#include "include/transformed_parameters_inner_deviation_exact_end.stan"
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
generated quantities {
#include "include/generated_quantities.stan"
}
