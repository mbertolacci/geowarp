int P_mean_total = P_mean_fixed + P_mean_random;
matrix[N, P_mean_total] X_mean = append_col(X_mean_fixed, X_mean_random);
vector[P_mean_total] alpha_Q_mu = rep_vector(0, P_mean_total);
vector[N] y_tilde = y - X_mean_fixed * alpha_mean;
