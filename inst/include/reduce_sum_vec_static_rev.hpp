#ifndef REDUCE_SUM_VEC_STATIC_REV_HPP
#define REDUCE_SUM_VEC_STATIC_REV_HPP

namespace geowarp {

namespace internal {

namespace {

template <
    typename ReduceFunction,
    typename... Args
>
class reduce_sum_vec_static_rev_vari : public stan::math::vari {
private:
    using args_tuple_t = std::tuple<Args...>;
    using args_values_tuple_t = std::tuple<
        stan::plain_type_t<
            decltype(value_of(std::declval<const Args&>()))
        >...
    >;
    using args_tuple_copy_t
        = std::tuple<decltype(deep_copy_vars(std::declval<Args>()))...>;
    using return_type_t = Eigen::Matrix<
        stan::math::var,
        Eigen::Dynamic,
        1
    >;

    struct recursive_reducer_forwards {
        std::stringstream msgs_;

        // Inputs
        const args_tuple_copy_t * const args_tuple_;

        // Outputs
        Eigen::VectorXd sum_;
        std::map<int, std::shared_ptr<ScopedChainableStack>> stacks_;
        std::map<int, std::shared_ptr<args_tuple_copy_t>> args_tuples_;
        std::map<int, std::shared_ptr<return_type_t>> outputs_;

        recursive_reducer_forwards(
            int output_rows,
            const args_tuple_copy_t * const args_tuple
        ) : args_tuple_(args_tuple),
            sum_(Eigen::VectorXd::Zero(output_rows)) {}

        recursive_reducer_forwards(
            recursive_reducer_forwards& other,
            tbb::split
        ) : args_tuple_(other.args_tuple_),
            sum_(Eigen::VectorXd::Zero(other.sum_.size())) {}

        inline void operator()(const tbb::blocked_range<size_t>& r) {
            if (r.empty()) {
                return;
            }

            int start = r.begin();
            stacks_[start] = std::make_shared<ScopedChainableStack>();

            (*stacks_[start]).execute([&]() {
                stan::math::apply([&](auto&&... args) {
                    args_tuples_[start] = std::make_shared<args_tuple_copy_t>(
                        deep_copy_vars(args)...
                    );
                }, *args_tuple_);

                stan::math::apply([&](auto && ... args) {
                    outputs_[start] = std::make_shared<return_type_t>(ReduceFunction()(
                        r.begin(),
                        r.end(),
                        &msgs_,
                        args...
                    ));
                }, *args_tuples_[start]);

                auto& output = *outputs_[start];
                for (int i = 0; i < sum_.rows(); ++i) {
                    sum_[i] += value_of(output[i]);
                }
            });
        }

        inline void join(const recursive_reducer_forwards& rhs) {
            sum_ += rhs.sum_;
            stacks_.insert(rhs.stacks_.begin(), rhs.stacks_.end());
            args_tuples_.insert(rhs.args_tuples_.begin(), rhs.args_tuples_.end());
            outputs_.insert(rhs.outputs_.begin(), rhs.outputs_.end());
            msgs_ << rhs.msgs_.str();
        }
    };

    struct recursive_reducer_backwards {
        std::stringstream msgs_;

        // Inputs
        std::map<int, std::shared_ptr<ScopedChainableStack>>& stacks_;
        std::map<int, std::shared_ptr<args_tuple_copy_t>>& args_tuples_;
        std::map<int, std::shared_ptr<return_type_t>>& outputs_;
        const Eigen::VectorXd& output_adjoints_;

        // Outputs
        Eigen::VectorXd args_adjoints_;

        recursive_reducer_backwards(
            int output_rows,
            int num_vars_args,
            std::map<int, std::shared_ptr<ScopedChainableStack>>& stacks,
            std::map<int, std::shared_ptr<args_tuple_copy_t>>& args_tuples,
            std::map<int, std::shared_ptr<return_type_t>>& outputs,
            const Eigen::VectorXd& output_adjoints
        ) : stacks_(stacks),
            args_tuples_(args_tuples),
            outputs_(outputs),
            output_adjoints_(output_adjoints),
            args_adjoints_(Eigen::VectorXd::Zero(num_vars_args)) {}

        recursive_reducer_backwards(
            recursive_reducer_backwards& other,
            tbb::split
        ) : stacks_(other.stacks_),
            args_tuples_(other.args_tuples_),
            outputs_(other.outputs_),
            output_adjoints_(other.output_adjoints_),
            args_adjoints_(Eigen::VectorXd::Zero(other.args_adjoints_.size())) {}

        inline void operator()(const tbb::blocked_range<size_t>& r) {
            if (r.empty()) {
                return;
            }

            auto& current_stack = *stacks_[r.begin()];
            auto& current_args_tuple = *args_tuples_[r.begin()];
            auto& current_output = *outputs_[r.begin()];

            current_stack.execute([&]() {
                for (int i = 0; i < current_output.size(); ++i) {
                    current_output[i].vi_->adj_ = output_adjoints_[i];
                }
                grad();
            });

            stan::math::apply(
                [&](auto&&... args) {
                    accumulate_adjoints(args_adjoints_.data(), args...);
                },
                current_args_tuple
            );
        }

        inline void join(const recursive_reducer_backwards& rhs) {
            args_adjoints_ += rhs.args_adjoints_;
            msgs_ << rhs.msgs_.str();
        }
    };

    struct reduce_sum_vec_static_rev_alloc : public chainable_alloc {
        virtual ~reduce_sum_vec_static_rev_alloc() {}

        args_tuple_copy_t args_tuple_copy_;
        recursive_reducer_forwards forwards_worker_;

        reduce_sum_vec_static_rev_alloc(
            int output_rows,
            Args&&... args
        ) : args_tuple_copy_(deep_copy_vars(args)...),
            forwards_worker_(output_rows, &args_tuple_copy_) {}
    };

    const int start_;
    const int end_;
    const int output_rows_;
    const int grainsize_;
    std::ostream *msgs_;

    const int num_vars_args_;
    vari **args_varis_;
    reduce_sum_vec_static_rev_alloc *alloc_;

public:
    stan::math::vari **output_;

    reduce_sum_vec_static_rev_vari(
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
        args_varis_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            num_vars_args_
        )),
        alloc_(new reduce_sum_vec_static_rev_alloc(output_rows, args...)),
        output_(
            ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari * >(output_rows_)
        ) {
        save_varis(args_varis_, args...);

        tbb::static_partitioner partitioner;
        tbb::parallel_deterministic_reduce(
            tbb::blocked_range<std::size_t>(start_, end_, grainsize_),
            alloc_->forwards_worker_,
            partitioner
        );

        if (msgs_) {
          (*msgs_) << alloc_->forwards_worker_.msgs_.str();
        }
        for (size_t i = 0; i < output_rows_; ++i) {
            output_[i] = new stan::math::vari(alloc_->forwards_worker_.sum_[i], false);
        }
    }

    virtual void chain() {
        Eigen::VectorXd output_adjoints(output_rows_);
        for (int i = 0; i < output_rows_; ++i) {
            output_adjoints[i] = output_[i]->adj_;
        }

        recursive_reducer_backwards backwards_worker(
            output_rows_,
            num_vars_args_,
            alloc_->forwards_worker_.stacks_,
            alloc_->forwards_worker_.args_tuples_,
            alloc_->forwards_worker_.outputs_,
            output_adjoints
        );

        tbb::static_partitioner partitioner;
        tbb::parallel_deterministic_reduce(
            tbb::blocked_range<std::size_t>(start_, end_, grainsize_),
            backwards_worker,
            partitioner
        );

        if (msgs_) {
          (*msgs_) << backwards_worker.msgs_.str();
        }
        for (int i = 0; i < num_vars_args_; ++i) {
            args_varis_[i]->adj_ += backwards_worker.args_adjoints_[i];
        }
    }
};

}

template <
    typename ReduceFunction,
    typename... Args
> struct reduce_sum_vec_impl<
    ReduceFunction,
    stan::math::var,
    false,
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
        auto vari = new reduce_sum_vec_static_rev_vari<
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

#endif  // REDUCE_SUM_VEC_STATIC_REV_HPP
