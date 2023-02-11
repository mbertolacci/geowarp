N_current_block = block_last_index[i] - current_block_start + 1;
N_parents_current_block = N_current_block - block_N_responses[i];

indices_current_block[:N_current_block] = block_indices[
  current_block_start:block_last_index[i]
];
