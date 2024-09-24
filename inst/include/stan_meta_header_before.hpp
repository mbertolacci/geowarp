#ifndef STAN_META_HEADER_BEFORE_HPP
#define STAN_META_HEADER_BEFORE_HPP

template <typename T_ss, typename T_sd, typename T_x>
Eigen::Matrix<
  stan::promote_args_t<T_ss, stan::value_type_t<T_sd>, T_x>,
  Eigen::Dynamic,
  Eigen::Dynamic
>
geowarp_process_covariance(
  const T_ss& sigma_squared_nugget,
  const T_sd& deviation_sd,
  const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1> >& x,
  const double& smoothness,
  std::ostream* pstream__
);

template <typename T_ss, typename T_sd, typename T_x>
Eigen::Matrix<
  stan::promote_args_t<T_ss, stan::value_type_t<T_sd>, T_x>,
  Eigen::Dynamic,
  Eigen::Dynamic
>
geowarp_process_covariance_1d(
    const T_ss& sigma_squared_nugget,
    const T_sd& deviation_sd_arg__,
    const std::vector<T_x>& x,
    const double& smoothness,
    std::ostream* pstream__
);

template <
    typename T3__,
    typename T4__,
    typename T5__,
    typename T6__,
    typename T7__,
    typename T8__
>
Eigen::Matrix<
    stan::promote_args_t<
        T3__,
        T4__,
        stan::value_type_t<T5__>,
        stan::value_type_t<T6__>,
        stan::value_type_t<T7__>,
        stan::promote_args_t<stan::value_type_t<T8__>>
    >,
    -1,
    1
>
geowarp_vecchia_reduce_sum_full(
    const int& start,
    const int& end,
    const int& grainsize,
    const std::vector<Eigen::Matrix<T3__, -1, 1>>& x_warped,
    const T4__& sigma_squared_nugget,
    const T5__& deviation_sd_arg__,
    const T6__& y_arg__,
    const T7__& y_tilde_arg__,
    const T8__& X_mean_arg__,
    const std::vector<int>& block_indices,
    const std::vector<int>& block_last_index,
    const std::vector<int>& block_N_responses,
    const std::vector<std::vector<int>>& indices_X_mean_non_zero,
    const std::vector<int>& P_X_mean_non_zero,
    const int& N_block_max,
    const double& smoothness,
    std::ostream* pstream__
);

template <
    typename T3__,
    typename T4__,
    typename T5__,
    typename T6__,
    typename T7__,
    typename T8__
>
Eigen::Matrix<
    stan::promote_args_t<
        T3__,
        T4__,
        stan::value_type_t<T5__>,
        stan::value_type_t<T6__>,
        stan::value_type_t<T7__>,
        stan::promote_args_t<stan::value_type_t<T8__>>
    >,
    -1,
    1
>
geowarp_vecchia_reduce_sum_vertical_only(
    const int& start,
    const int& end,
    const int& grainsize,
    const std::vector<T3__>& x_vertical_warped,
    const T4__& sigma_squared_nugget,
    const T5__& deviation_sd_arg__,
    const T6__& y_arg__,
    const T7__& y_tilde_arg__,
    const T8__& X_mean_arg__,
    const std::vector<int>& block_indices,
    const std::vector<int>& block_last_index,
    const std::vector<int>& block_N_responses,
    const std::vector<std::vector<int>>& indices_X_mean_non_zero,
    const std::vector<int>& P_X_mean_non_zero,
    const int& N_block_max,
    const double& smoothness,
    std::ostream* pstream__
);

#endif  // STAN_META_HEADER_BEFORE_HPP
