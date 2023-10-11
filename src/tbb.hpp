#ifndef TBB_HPP

int set_tbb_threads(int n_threads);
int get_tbb_threads();

class SetTbbThreads {
public:
    SetTbbThreads(int n_threads) {
        set_tbb_threads(n_threads);
    }

    ~SetTbbThreads() {
        set_tbb_threads(0);
    }
};

#endif
