#include "include/transformed_data_declarations.stan"

int N_block_max = get_N_block_max(block_last_index);
int indices_X_mean_non_zero[P_mean_total, N_blocks];
int P_X_mean_non_zero[N_blocks];
int current_block_start_a = 1;

#include "include/transformed_data_statements.stan"

for (i in 1:N_blocks) {
  int N_current_block = block_last_index[i] - current_block_start_a + 1;
  int indices_current_block[N_block_max];

  indices_current_block[:N_current_block] = block_indices[
    current_block_start_a:block_last_index[i]
  ];
  P_X_mean_non_zero[i] = 0;
  for (k in 1:P_mean_total) {
    if (max(fabs(X_mean[indices_current_block[:N_current_block], k])) != 0) {
      P_X_mean_non_zero[i] += 1;
      indices_X_mean_non_zero[P_X_mean_non_zero[i], i] = k;
    }
  }

  current_block_start_a = block_last_index[i] + 1;
}
