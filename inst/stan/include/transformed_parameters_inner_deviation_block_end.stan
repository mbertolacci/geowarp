L_block[:N_current_block, :N_current_block] = cholesky_decompose(
  K_block[:N_current_block, :N_current_block]
);
L_inv_X_block[:N_current_block, :P_X_mean_non_zero[i]] = mdivide_left_tri_low(
  L_block[:N_current_block, :N_current_block],
  X_mean[
    indices_current_block[:N_current_block],
    indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]
  ]
);
L_inv_y_block[:N_current_block] = mdivide_left_tri_low(
  L_block[:N_current_block, :N_current_block],
  y[indices_current_block[:N_current_block]]
);
L_inv_y_tilde_block[:N_current_block] = mdivide_left_tri_low(
  L_block[:N_current_block, :N_current_block],
  y_tilde[indices_current_block[:N_current_block]]
);

log_det += (
  2 * sum(log(diagonal(L_block[
    (N_parents_current_block+1):N_current_block,
    (N_parents_current_block+1):N_current_block
  ])))
);
y_tildet_Q_y_tilde += (
  sum(square(L_inv_y_tilde_block[(N_parents_current_block+1):N_current_block]))
);
Xt_Q_y[indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]] += (
  L_inv_X_block[(N_parents_current_block+1):N_current_block, :P_X_mean_non_zero[i]]'
  * L_inv_y_block[(N_parents_current_block+1):N_current_block]
);
Xt_Q_y_tilde[indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]] += (
  L_inv_X_block[(N_parents_current_block+1):N_current_block, :P_X_mean_non_zero[i]]'
  * L_inv_y_tilde_block[(N_parents_current_block+1):N_current_block]
);
Xt_Q_X[
  indices_X_mean_non_zero[:P_X_mean_non_zero[i], i],
  indices_X_mean_non_zero[:P_X_mean_non_zero[i], i]
] += (
  crossprod(L_inv_X_block[
    (N_parents_current_block+1):N_current_block,
    :P_X_mean_non_zero[i]
  ])
);

current_block_start = block_last_index[i] + 1;
