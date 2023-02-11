#include "include/transformed_data_declarations.stan"

int N_block_max = get_N_block_max(block_last_index);
int indices_X_mean_non_zero[P_mean_total, N_blocks];
matrix[N_block_max, P_mean_total] X_mean_non_zero[N_blocks];
int P_X_mean_non_zero[N_blocks];
int current_block_start_a = 1;

#include "include/transformed_data_statements.stan"

for (i in 1:N_blocks) {
  vector[N_block_max] X_mean_block_k;
  int N_current_block = block_last_index[i] - current_block_start_a + 1;
  int indices_current_block[N_block_max];

  indices_current_block[:N_current_block] = block_indices[
    current_block_start_a:block_last_index[i]
  ];
  P_X_mean_non_zero[i] = 0;
  for (k in 1:P_mean_total) {
    X_mean_block_k[:N_current_block] = X_mean[indices_current_block[:N_current_block], k];
    if (max(fabs(X_mean_block_k[:N_current_block])) != 0) {
      P_X_mean_non_zero[i] += 1;
      indices_X_mean_non_zero[P_X_mean_non_zero[i], i] = k;
      X_mean_non_zero[
        i,
        :N_current_block,
        P_X_mean_non_zero[i]
      ] = X_mean_block_k[:N_current_block];
    }
  }

  current_block_start_a = block_last_index[i] + 1;
}
