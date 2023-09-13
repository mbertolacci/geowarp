target += log_marginal;

if (tau_squared_mean_random_b > 0) {
  tau_squared_mean_random ~ inv_gamma(tau_squared_mean_random_a, tau_squared_mean_random_b);
} else {
  target += -log(tau_squared_mean_random);
}

if (sigma_squared_nugget_b > 0) {
  sigma_squared_nugget ~ inv_gamma(sigma_squared_nugget_a, sigma_squared_nugget_b);
} else {
  target += -log(sigma_squared_nugget);
}
