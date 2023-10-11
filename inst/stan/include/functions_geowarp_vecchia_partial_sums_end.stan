 L_block = cholesky_decompose(K_block);
  L_inv_X_block = mdivide_left_tri_low(
    L_block,
    X_mean[
      indices_current_block,
      indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]
    ]
  );
  L_inv_y_block = mdivide_left_tri_low(
    L_block,
    y[indices_current_block]
  );
  L_inv_y_tilde_block = mdivide_left_tri_low(
    L_block,
    y_tilde[indices_current_block]
  );
  log_det += (
    2 * sum(log(diagonal(L_block[
      (N_parents_current_block+1):,
      (N_parents_current_block+1):
    ])))
  );
  y_tildet_Q_y_tilde += (
    sum(square(L_inv_y_tilde_block[(N_parents_current_block + 1):]))
  );
  Xt_Q_y[indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]] += (
    L_inv_X_block[(N_parents_current_block + 1):, :P_X_mean_non_zero[i]]'
    * L_inv_y_block[(N_parents_current_block + 1):]
  );
  Xt_Q_y_tilde[indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]] += (
    L_inv_X_block[(N_parents_current_block + 1):, :P_X_mean_non_zero[i]]'
    * L_inv_y_tilde_block[(N_parents_current_block + 1):]
  );
  Xt_Q_X[
    indices_X_mean_non_zero[:P_X_mean_non_zero[i], i],
    indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]
  ] += (
    crossprod(L_inv_X_block[(N_parents_current_block + 1):, :])
  );
}

vector[2 + P_mean * P_mean + 2 * P_mean] output;
output[1] = log_det;
output[2] = y_tildet_Q_y_tilde;
output[3:(2 + P_mean * P_mean)] = to_vector(Xt_Q_X);
output[(3 + P_mean * P_mean):(2 + P_mean * P_mean + P_mean)] = to_vector(Xt_Q_y);
output[(3 + P_mean * P_mean + P_mean):] = to_vector(Xt_Q_y_tilde);
return output;
