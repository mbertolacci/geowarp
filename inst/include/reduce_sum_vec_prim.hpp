#ifndef REDUCE_SUM_VEC_PRIM_HPP
#define REDUCE_SUM_VEC_PRIM_HPP

namespace geowarp {

namespace internal {

// Specialisation for double type
template <
    typename ReduceFunction,
    bool Dynamic,
    typename... Args
> struct reduce_sum_vec_impl<
    ReduceFunction,
    double,
    Dynamic,
    Args...
> {
    using args_tuple_t = std::tuple<Args...>;

    struct recursive_reducer {
        args_tuple_t args_tuple_;

        std::stringstream msgs_;
        Eigen::VectorXd sum_;

        recursive_reducer(
            int output_rows,
            Args&&... args
        ) : args_tuple_(std::forward<Args>(args)...),
            sum_(Eigen::VectorXd::Zero(output_rows)) {}

        recursive_reducer(
            recursive_reducer& other,
            tbb::split
        ) : args_tuple_(other.args_tuple_),
            sum_(Eigen::VectorXd::Zero(other.sum_.size())) {}

        inline void operator()(const tbb::blocked_range<size_t>& r) {
            if (r.empty()) {
                return;
            }

            sum_ += stan::math::apply(
                [&](auto && ... args) {
                    return ReduceFunction()(
                        r.begin(),
                        r.end(),
                        &msgs_,
                        args...
                    );
                },
                args_tuple_
            );
        }

        inline void join(const recursive_reducer& rhs) {
            sum_ += rhs.sum_;
            msgs_ << rhs.msgs_.str();
        }
    };

    Eigen::VectorXd operator()(
        int start,
        int end,
        int output_rows,
        int grainsize,
        std::ostream* msgs,
        Args&&... args
    ) const {
        recursive_reducer worker(output_rows, std::forward<Args>(args)...);

        if (Dynamic) {
            tbb::parallel_reduce(
                tbb::blocked_range<std::size_t>(start, end, grainsize),
                worker
            );
        } else {
            tbb::static_partitioner partitioner;
            tbb::parallel_deterministic_reduce(
                tbb::blocked_range<std::size_t>(start, end, grainsize),
                worker,
                partitioner
            );
        }

        if (msgs) {
            (*msgs) << worker.msgs_.str();
        }

        return worker.sum_;
    }
};

}  // namespace internal

}  // namespace geowarp

#endif  // REDUCE_SUM_VEC_PRIM_HPP
