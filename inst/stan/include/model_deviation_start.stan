eta_deviation ~ multi_normal_prec(eta_deviation_mean, eta_deviation_precision);

if (P_deviation_random > 0) {
  zeta_deviation ~ multi_normal_prec(
    rep_vector(0, P_deviation_random),
    exp1d_precision(
      P_deviation_random,
      delta_deviation_random,
      ell_deviation_random,
      tau_squared_deviation_random
    )
  );
}

if (ell_deviation_random_scale > 0) {
  ell_deviation_random ~ normal(0, ell_deviation_random_scale);
}

if (tau_squared_deviation_random_b > 0) {
  tau_squared_deviation_random ~ inv_gamma(tau_squared_deviation_random_a, tau_squared_deviation_random_b);
} else {
  target += -log(tau_squared_deviation_random);
}
