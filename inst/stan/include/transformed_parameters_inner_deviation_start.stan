vector[N] deviation_log_variance;
vector[N] deviation_sd;
matrix[N_block_max, N_block_max] K_block;
matrix[N_block_max, N_block_max] L_block;
matrix[N_block_max, P_mean_total] L_inv_X_block;
vector[N_block_max] L_inv_y_block;
vector[N_block_max] L_inv_y_tilde_block;
int current_block_start = 1;
int N_current_block;
int N_parents_current_block;
int indices_current_block[N_block_max];

deviation_log_variance = X_deviation_fixed * eta_deviation;
if (P_deviation_random > 0) {
  deviation_log_variance += X_deviation_random * zeta_deviation;
}
deviation_sd = exp(0.5 * deviation_log_variance);
