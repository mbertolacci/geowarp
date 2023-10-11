#ifndef STAN_META_HEADER_HPP
#define STAN_META_HEADER_HPP

#include <geowarp_process_covariance.hpp>
#include <reduce_sum_vec.hpp>
#include <reduce_sum_vec_prim.hpp>
#include <reduce_sum_vec_dynamic_rev.hpp>
#include <reduce_sum_vec_static_rev.hpp>

struct geowarp_vecchia_partial_sums_full_rsvfunctor {
    template <
        typename T2__,
        typename T3__,
        typename T4__,
        typename T5__,
        typename T6__,
        typename T7__
    >
    Eigen::Matrix<
        stan::promote_args_t<
            T2__,
            T3__,
            stan::value_type_t<T4__>,
            stan::value_type_t<T5__>,
            stan::value_type_t<T6__>,
            stan::promote_args_t<stan::value_type_t<T7__>>
        >,
        -1,
        1
    >
    operator() (
        const int& start,
        const int& end,
        std::ostream* pstream__,
        const std::vector<Eigen::Matrix<T2__, -1, 1>>& x_warped,
        const T3__& sigma_squared_nugget,
        const T4__& deviation_sd_arg__,
        const T5__& y_arg__,
        const T6__& y_tilde_arg__,
        const T7__& X_mean_arg__,
        const std::vector<int>& block_indices,
        const std::vector<int>& block_last_index,
        const std::vector<int>& block_N_responses,
        const std::vector<std::vector<int>>& indices_X_mean_non_zero,
        const std::vector<int>& P_X_mean_non_zero,
        const int& N_block_max,
        const double& smoothness
    ) {
        return geowarp_vecchia_partial_sums_full(
            start + 1,
            end,
            x_warped,
            sigma_squared_nugget,
            deviation_sd_arg__,
            y_arg__,
            y_tilde_arg__,
            X_mean_arg__,
            block_indices,
            block_last_index,
            block_N_responses,
            indices_X_mean_non_zero,
            P_X_mean_non_zero,
            N_block_max,
            smoothness,
            pstream__
        );
    }
};

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
) {
    int P = X_mean_arg__.cols();
    return geowarp::reduce_sum_vec_dynamic<geowarp_vecchia_partial_sums_full_rsvfunctor>(
        start,
        end,
        2 + P * P + 2 * P,
        grainsize,
        pstream__,
        x_warped,
        sigma_squared_nugget,
        deviation_sd_arg__,
        y_arg__,
        y_tilde_arg__,
        X_mean_arg__,
        block_indices,
        block_last_index,
        block_N_responses,
        indices_X_mean_non_zero,
        P_X_mean_non_zero,
        N_block_max,
        smoothness
    );
}

struct geowarp_vecchia_partial_sums_vertical_only_rsvfunctor {
    template <
        typename T2__,
        typename T3__,
        typename T4__,
        typename T5__,
        typename T6__,
        typename T7__
    >
    Eigen::Matrix<
        stan::promote_args_t<
            T2__,
            T3__,
            stan::value_type_t<T4__>,
            stan::value_type_t<T5__>,
            stan::value_type_t<T6__>,
            stan::promote_args_t<stan::value_type_t<T7__>>
        >,
        -1,
        1
    >
    operator() (
        const int& start,
        const int& end,
        std::ostream* pstream__,
        const std::vector<T2__>& x_vertical_warped,
        const T3__& sigma_squared_nugget,
        const T4__& deviation_sd_arg__,
        const T5__& y_arg__,
        const T6__& y_tilde_arg__,
        const T7__& X_mean_arg__,
        const std::vector<int>& block_indices,
        const std::vector<int>& block_last_index,
        const std::vector<int>& block_N_responses,
        const std::vector<std::vector<int>>& indices_X_mean_non_zero,
        const std::vector<int>& P_X_mean_non_zero,
        const int& N_block_max,
        const double& smoothness
    ) {
        return geowarp_vecchia_partial_sums_vertical_only(
            start + 1,
            end,
            x_vertical_warped,
            sigma_squared_nugget,
            deviation_sd_arg__,
            y_arg__,
            y_tilde_arg__,
            X_mean_arg__,
            block_indices,
            block_last_index,
            block_N_responses,
            indices_X_mean_non_zero,
            P_X_mean_non_zero,
            N_block_max,
            smoothness,
            pstream__
        );
    }
};

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
) {
    int P = X_mean_arg__.cols();
    return geowarp::reduce_sum_vec_dynamic<geowarp_vecchia_partial_sums_vertical_only_rsvfunctor>(
        start,
        end,
        2 + P * P + 2 * P,
        grainsize,
        pstream__,
        x_vertical_warped,
        sigma_squared_nugget,
        deviation_sd_arg__,
        y_arg__,
        y_tilde_arg__,
        X_mean_arg__,
        block_indices,
        block_last_index,
        block_N_responses,
        indices_X_mean_non_zero,
        P_X_mean_non_zero,
        N_block_max,
        smoothness
    );
}

#endif  // STAN_META_HEADER_HPP
