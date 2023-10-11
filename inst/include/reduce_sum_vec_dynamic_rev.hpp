#ifndef REDUCE_SUM_VEC_DYNAMIC_REV_HPP
#define REDUCE_SUM_VEC_DYNAMIC_REV_HPP

namespace geowarp {

namespace internal {

namespace {

template <
    typename ReduceFunction,
    typename... Args
>
class reduce_sum_vec_dynamic_rev_vari : public stan::math::vari {
private:
    using args_tuple_t = std::tuple<Args...>;
    using args_tuple_copy_t = std::tuple<
        decltype(deep_copy_vars(std::declval<Args>()))...
    >;
    using args_values_tuple_t = std::tuple<
        stan::plain_type_t<
            decltype(value_of(std::declval<const Args&>()))
        >...
    >;

    // This reducer stores the arguments as raw values, which will call the
    // specialisation of ReduceFunction for double types. This is usually
    // faster. It accumulates the sum of the outputs.
    struct recursive_reducer_forwards {
        const args_values_tuple_t * const args_values_tuple_;

        std::stringstream msgs_;
        Eigen::VectorXd sum_;

        recursive_reducer_forwards(
            int output_rows,
            const args_values_tuple_t * const args_values_tuple
        ) : args_values_tuple_(args_values_tuple),
            sum_(Eigen::VectorXd::Zero(output_rows)) {}

        recursive_reducer_forwards(
            recursive_reducer_forwards& other,
            tbb::split
        ) : args_values_tuple_(other.args_values_tuple_),
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
                *args_values_tuple_
            );
        }

        inline void join(const recursive_reducer_forwards& rhs) {
            sum_ += rhs.sum_;
            msgs_ << rhs.msgs_.str();
        }
    };

    // This reducer uses the original inputs, and accumulates the adjoints of
    // the inputs
    struct recursive_reducer_backwards {
        std::stringstream msgs_;

        // Inputs
        const args_tuple_copy_t * const args_tuple_;
        const Eigen::VectorXd& output_adjoints_;

        // Internal storage
        std::unique_ptr<ScopedChainableStack> stack_;
        std::unique_ptr<args_tuple_copy_t> local_args_tuple_copy_;

        // Outputs
        Eigen::VectorXd args_adjoints_;

        recursive_reducer_backwards(
            int output_rows,
            int num_vars_args,
            const Eigen::VectorXd& output_adjoints,
            const args_tuple_copy_t * const args_tuple
        ) : args_tuple_(args_tuple),
            output_adjoints_(output_adjoints),
            stack_(nullptr),
            local_args_tuple_copy_(nullptr),
            args_adjoints_(Eigen::VectorXd::Zero(num_vars_args)) {}

        recursive_reducer_backwards(
            recursive_reducer_backwards& other,
            tbb::split
        ) : args_tuple_(other.args_tuple_),
            output_adjoints_(other.output_adjoints_),
            stack_(nullptr),
            local_args_tuple_copy_(nullptr),
            args_adjoints_(Eigen::VectorXd::Zero(other.args_adjoints_.size())) {}

        inline void operator()(const tbb::blocked_range<size_t>& r) {
            if (r.empty()) {
                return;
            }

            if (!stack_) {
                stack_ = std::make_unique<ScopedChainableStack>();
            }

            stack_->execute([&]() {
                if (!local_args_tuple_copy_) {
                    local_args_tuple_copy_ = stan::math::apply(
                        [&](auto&&... args) {
                            return std::make_unique<args_tuple_copy_t>(
                                deep_copy_vars(args)...
                            );
                        },
                        *args_tuple_
                    );
                } else {
                    set_zero_all_adjoints();
                }
            });

            const nested_rev_autodiff begin_nest;

            auto& args_tuple_local = *local_args_tuple_copy_;
            auto output = stan::math::apply(
                [&](auto && ... args) {
                    return ReduceFunction()(
                        r.begin(),
                        r.end(),
                        &msgs_,
                        args...
                    );
                },
                args_tuple_local
            );

            for (int i = 0; i < output.rows(); ++i) {
                output[i].vi_->adj_ = output_adjoints_[i];
            }
            grad();

            stan::math::apply(
                [&](auto&&... args) {
                    accumulate_adjoints(args_adjoints_.data(), args...);
                },
                args_tuple_local
            );
        }

        inline void join(const recursive_reducer_backwards& rhs) {
            args_adjoints_ += rhs.args_adjoints_;
            msgs_ << rhs.msgs_.str();
        }
    };

    struct reduce_sum_vec_dynamic_rev_alloc : public chainable_alloc {
        virtual ~reduce_sum_vec_dynamic_rev_alloc() {}

        args_tuple_copy_t args_tuple_copy_;

        reduce_sum_vec_dynamic_rev_alloc(
            Args&&... args
        ) : args_tuple_copy_(deep_copy_vars(args)...) {}
    };

    const int start_;
    const int end_;
    const int output_rows_;
    const int grainsize_;
    std::ostream* msgs_;

    const int num_vars_args_;
    vari **args_varis_;
    reduce_sum_vec_dynamic_rev_alloc *alloc_;

public:
    stan::math::vari **output_;

    reduce_sum_vec_dynamic_rev_vari(
        int start,
        int end,
        int output_rows,
        int grainsize,
        std::ostream* msgs,
        Args&&... args
    ) : vari(0.0),
        start_(start),
        end_(end),
        output_rows_(output_rows),
        grainsize_(grainsize),
        msgs_(msgs),
        num_vars_args_(count_vars(std::forward<Args>(args)...)),
        args_varis_(ChainableStack::instance_->memalloc_.alloc_array<vari *>(
            num_vars_args_
        )),
        alloc_(new reduce_sum_vec_dynamic_rev_alloc(args...)),
        output_(
            ChainableStack::instance_->memalloc_.alloc_array<vari *>(output_rows_)
        ) {
        save_varis(args_varis_, args...);

        std::unique_ptr<args_values_tuple_t> args_value_tuple = std::make_unique<args_values_tuple_t>(
            // TODO: I think this currently copies arguments that are
            // already values; it would be nice/faster to save their
            // references
            value_of(std::forward<Args>(args))...
        );
        recursive_reducer_forwards worker(
            output_rows_,
            args_value_tuple.get()
        );

        tbb::parallel_reduce(
            tbb::blocked_range<std::size_t>(start_, end_, grainsize_),
            worker
        );

        if (msgs_) {
          (*msgs_) << worker.msgs_.str();
        }
        for (size_t i = 0; i < output_rows_; ++i) {
            output_[i] = new stan::math::vari(worker.sum_[i], false);
        }
    }

    virtual void chain() {
        Eigen::VectorXd output_adjoints(output_rows_);
        for (int i = 0; i < output_rows_; ++i) {
            output_adjoints[i] = output_[i]->adj_;
        }

        recursive_reducer_backwards worker(
            output_rows_,
            num_vars_args_,
            output_adjoints,
            &alloc_->args_tuple_copy_
        );

        tbb::parallel_reduce(
            tbb::blocked_range<std::size_t>(start_, end_, grainsize_),
            worker
        );

        if (msgs_) {
          (*msgs_) << worker.msgs_.str();
        }
        for (int i = 0; i < num_vars_args_; ++i) {
            args_varis_[i]->adj_ += worker.args_adjoints_[i];
        }
    }
};

}  // empty namespace

template <
    typename ReduceFunction,
    typename... Args
> struct reduce_sum_vec_impl<
    ReduceFunction,
    stan::math::var,
    true,
    Args...
> {
    using ReturnType = Eigen::Matrix<
        stan::math::var,
        Eigen::Dynamic,
        1
    >;

    ReturnType operator()(
        int start,
        int end,
        int output_rows,
        int grainsize,
        std::ostream* msgs,
        Args&&... args
    ) const {
        auto vari = new reduce_sum_vec_dynamic_rev_vari<
            ReduceFunction,
            stan::ref_type_t<Args&&>...
        >(
            start,
            end,
            output_rows,
            grainsize,
            msgs,
            std::forward<Args>(args)...
        );

        ReturnType output(output_rows);
        for (size_t i = 0; i < output_rows; ++i) {
            output[i].vi_ = vari->output_[i];
        }
        return output;
    }
};

}  // namespace internal

}  // namespace geowarp

#endif  // REDUCE_SUM_VEC_DYNAMIC_REV_HPP
