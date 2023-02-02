vector[N] deviation_log_variance;
vector[N] deviation_sd;
matrix[N, N] L;
matrix[N, P_mean_total] L_inv_X;
vector[N] L_inv_y;
vector[N] L_inv_y_tilde;

deviation_log_variance = X_deviation_fixed * eta_deviation;
if (P_deviation_random > 0) {
  deviation_log_variance += X_deviation_random * zeta_deviation;
}
deviation_sd = exp(0.5 * deviation_log_variance);
