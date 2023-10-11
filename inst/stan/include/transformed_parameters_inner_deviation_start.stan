vector[N] deviation_log_variance;
vector[N] deviation_sd;
vector[
  2 + P_mean_total * P_mean_total + 2 * P_mean_total
] components;

deviation_log_variance = X_deviation_fixed * eta_deviation;
if (P_deviation_random > 0) {
  deviation_log_variance += X_deviation_random * zeta_deviation;
}
deviation_sd = exp(0.5 * deviation_log_variance);
