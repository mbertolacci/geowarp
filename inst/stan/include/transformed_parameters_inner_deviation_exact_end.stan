for (i in 1:N) {
  K[i, i] = K[i, i] + sigma_squared_nugget;
}
L = cholesky_decompose(K);

L_inv_X = mdivide_left_tri_low(L, X_mean);
L_inv_y = mdivide_left_tri_low(L, y);
L_inv_y_tilde = mdivide_left_tri_low(L, y_tilde);

log_det = 2 * sum(log(diagonal(L)));
yt_Q_y = sum(square(L_inv_y));
y_tildet_Q_y_tilde = sum(square(L_inv_y_tilde));
Xt_Q_X = crossprod(L_inv_X);
Xt_Q_y = L_inv_X' * L_inv_y;
Xt_Q_y_tilde = L_inv_X' * L_inv_y_tilde;
