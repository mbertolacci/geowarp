#include <tbb/global_control.h>
#include <tbb/task_arena.h>

#include "tbb.hpp"

std::unique_ptr<tbb::global_control> control;

// [[Rcpp::export(name = ".set_tbb_threads", rng = false)]]
int set_tbb_threads(int n_threads) {
    int old_n_threads = -1;
    if (control) {
        old_n_threads = control->active_value(tbb::global_control::max_allowed_parallelism);
    }
    control.reset();
    if (n_threads > 0) {
        control = std::make_unique<tbb::global_control>(
            tbb::global_control::max_allowed_parallelism, n_threads
        );
    }
    return old_n_threads;
}

// [[Rcpp::export(name = ".get_tbb_threads", rng = false)]]
int get_tbb_threads() {
    if (control) {
        return control->active_value(tbb::global_control::max_allowed_parallelism);
    }
    return tbb::this_task_arena::max_concurrency();
}
