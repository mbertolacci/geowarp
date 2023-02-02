eta_deviation ~ multi_normal_prec(eta_deviation_mean, eta_deviation_precision);

if (P_deviation_random > 0) {
  // matrix[
  //   P_deviation_random,
  //   P_deviation_random
  // ] zeta_deviation_precision_L = cholesky_decompose(exp1d_precision(
  //   P_deviation_random,
  //   1,
  //   ell_deviation_random,
  //   tau_squared_deviation_random_b
  // ));
  // target += sum(log(diagonal(zeta_deviation_precision_L)));
  // target += -0.5 * sum(square(
  //   zeta_deviation_precision_L' * zeta_deviation
  // ));

  zeta_deviation ~ multi_normal_prec(
    rep_vector(0, P_deviation_random),
    // rw1d_precision(P_deviation_random, tau_squared_deviation_random)
    exp1d_precision(
      P_deviation_random,
      delta_deviation_random,
      ell_deviation_random,
      tau_squared_deviation_random
    )
  );

  // zeta_deviation ~ multi_normal(
  //   rep_vector(0, P_deviation_random),
  //   exp1d_covariance(
  //     P_deviation_random,
  //     1,
  //     ell_deviation_random,
  //     tau_squared_deviation_random,
  //     tau_squared_deviation_random2
  //   )
  // );
}

if (ell_deviation_random_scale > 0) {
  ell_deviation_random ~ normal(0, ell_deviation_random_scale);
}

if (tau_squared_deviation_random_scale > 0) {
  tau_squared_deviation_random ~ normal(0, tau_squared_deviation_random_scale);
} else {
  target += -log(tau_squared_deviation_random);
}
