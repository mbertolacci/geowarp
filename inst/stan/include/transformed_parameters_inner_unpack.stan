log_det = components[1];
y_tildet_Q_y_tilde = components[2];
Xt_Q_X = to_matrix(
  components[3:(2 + P_mean_total * P_mean_total)],
  P_mean_total,
  P_mean_total
);
Xt_Q_y = components[
  (3 + P_mean_total * P_mean_total):(
    2 + P_mean_total * P_mean_total + P_mean_total
  )
];
Xt_Q_y_tilde = components[
  (3 + P_mean_total * P_mean_total + P_mean_total):
];
