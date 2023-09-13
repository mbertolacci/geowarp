alpha_beta_prior_precision[:P_mean_fixed, :P_mean_fixed] = alpha_precision;
if (P_mean_total > P_mean_fixed) {
  alpha_beta_prior_precision[
    (P_mean_fixed+1):P_mean_total,
    (P_mean_fixed+1):P_mean_total
  ] = rw1d_precision(
    P_mean_random,
    tau_squared_mean_random
  );
}

alpha_beta_precision = Xt_Q_X + alpha_beta_prior_precision;

L_alpha_beta_prior_precision = cholesky_decompose(alpha_beta_prior_precision);
L_alpha_beta_precision = cholesky_decompose(alpha_beta_precision);
alpha_beta_rhs = Xt_Q_y + alpha_Q_mu;
alpha_beta_hat = chol_solve_L_b(L_alpha_beta_precision, alpha_beta_rhs);

// This uses the Woodbury matrix identity
log_marginal = -0.5 * (
  2 * sum(log(diagonal(L_alpha_beta_precision)))
  - 2 * sum(log(diagonal(L_alpha_beta_prior_precision)))
  + log_det
  + y_tildet_Q_y_tilde
  - sum(square(mdivide_left_tri_low(
    L_alpha_beta_precision,
    Xt_Q_y_tilde
  )))
);
