int P_mean = cols(X_mean);

real log_det = 0;
real y_tildet_Q_y_tilde = 0;
matrix[P_mean, P_mean] Xt_Q_X = rep_matrix(0, P_mean, P_mean);
vector[P_mean] Xt_Q_y = rep_vector(0, P_mean);
vector[P_mean] Xt_Q_y_tilde = rep_vector(0, P_mean);

for (i in start:end) {
  int block_start = get_block_start(i, block_last_index);
  int N_current_block = block_last_index[i] - block_start + 1;
  matrix[N_current_block, N_current_block] K_block;
  matrix[N_current_block, N_current_block] L_block;
  matrix[N_current_block, P_X_mean_non_zero[i]] L_inv_X_block;
  vector[N_current_block] L_inv_y_block;
  vector[N_current_block] L_inv_y_tilde_block;
  int indices_current_block[N_current_block] = block_indices[
    block_start:block_last_index[i]
  ];
  int N_parents_current_block = N_current_block - block_N_responses[i];
