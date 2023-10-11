#ifndef REDUCE_SUM_VEC_HPP
#define REDUCE_SUM_VEC_HPP

namespace geowarp {

namespace internal {

template <
    typename ReduceFunction,
    typename BaseReturnType,
    bool Dynamic,
    typename... Args
>
struct reduce_sum_vec_impl;

}  // namespace internal

template <
    typename ReduceFunction,
    typename... Args
>
inline auto reduce_sum_vec_dynamic(
    int start,
    int end,
    int output_rows,
    int grainsize,
    std::ostream* msgs,
    Args&&... args
) {
    // Grab the return type of the ReduceFunction
    using ReturnType = decltype(
        ReduceFunction()(
            start,
            end,
            msgs,
            std::forward<Args>(args)...
        )
    );
    using BaseReturnType = typename ReturnType::Scalar;

    // TODO: add more checks
    check_positive("reduce_sum_vec_dynamic", "output_rows", output_rows);
    check_positive("reduce_sum_vec_dynamic", "grainsize", grainsize);

#ifdef STAN_THREADS
    return internal::reduce_sum_vec_impl<
        ReduceFunction,
        BaseReturnType,
        true,
        stan::ref_type_t<Args&&>...
    >()(
        start - 1,
        end,
        output_rows,
        grainsize,
        msgs,
        std::forward<Args>(args)...
    );
#else
    return ReduceFunction()(
        start - 1,
        end,
        msgs,
        std::forward<Args>(args)...
    );
#endif
}

template <
    typename ReduceFunction,
    typename... Args
>
inline auto reduce_sum_vec_static(
    int start,
    int end,
    int output_rows,
    int grainsize,
    std::ostream* msgs,
    Args&&... args
) {
    // Grab the return type of the ReduceFunction
    using ReturnType = decltype(
        ReduceFunction()(
            start,
            end,
            msgs,
            std::forward<Args>(args)...
        )
    );
    using BaseReturnType = typename ReturnType::Scalar;

    // TODO: add more checks
    check_positive("reduce_sum_vec_static", "output_rows", output_rows);
    check_positive("reduce_sum_vec_static", "grainsize", grainsize);

#ifdef STAN_THREADS
    return internal::reduce_sum_vec_impl<
        ReduceFunction,
        BaseReturnType,
        false,
        stan::ref_type_t<Args&&>...
    >()(
        start - 1,
        end,
        output_rows,
        grainsize,
        msgs,
        std::forward<Args>(args)...
    );
#else
    return ReduceFunction()(
        start - 1,
        end,
        msgs,
        std::forward<Args>(args)...
    );
#endif
}

}  // namespace geowarp

#endif  // REDUCE_SUM_VEC_HPP
