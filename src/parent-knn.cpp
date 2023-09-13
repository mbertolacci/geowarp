#include <Rcpp.h>
#include <RcppEigen.h>

#include "nanoflann.hpp"

template <
    typename _DistanceType, typename _IndexType = size_t,
    typename _CountType = size_t
>
class ParentKNNResultSet {
public:
    using DistanceType = _DistanceType;
    using IndexType    = _IndexType;
    using CountType    = _CountType;

private:
    IndexType target;
    std::vector<IndexType> indices;
    std::vector<DistanceType> dists;
    CountType capacity;
    CountType count;

public:
    explicit ParentKNNResultSet(CountType capacity_)
        : indices(capacity_),
          dists(capacity_),
          capacity(capacity_),
          count(0) { }

    void reset(IndexType target_) {
        target = target_;
        count = 0;
        dists[capacity - 1] = std::numeric_limits<DistanceType>::max();
    }

    CountType size() const {
        return count;
    }

    bool full() const {
        return count == capacity;
    }

    IndexType index(int i) {
        return indices[i];
    }

    /**
     * Called during search to add an element matching the criteria.
     * @return true if the search should be continued, false if the results are
     * sufficient
     */
    bool addPoint(DistanceType dist, IndexType index) {
        if (index > target) {
            // Continue the search, rejecting the point
            return true;
        }
        CountType i;
        for (i = count; i > 0; --i) {
            if (dists[i - 1] > dist) {
                if (i < capacity) {
                    dists[i] = dists[i - 1];
                    indices[i] = indices[i - 1];
                }
            } else {
                break;
            }
        }
        if (i < capacity) {
            dists[i]   = dist;
            indices[i] = index;
        }
        if (count < capacity) count++;

        // tell caller that the search shall continue
        return true;
    }

    DistanceType worstDist() const {
        return dists[capacity - 1];
    }
};

// [[Rcpp::export(name = ".parent_knn", rng = false)]]
Rcpp::IntegerMatrix parent_knn(
    const Eigen::MatrixXd& x,
    int n_parents,
    int leaf_size = 40
) {
    typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> kd_tree_type;

    int k = n_parents + 1;
    kd_tree_type kdTree(x.cols(), x, leaf_size);
    Rcpp::IntegerMatrix output(x.rows(), k);
    output.fill(Rcpp::IntegerVector::get_na());

    ParentKNNResultSet<double> resultSet(k);

    for (int i = 0; i < x.rows(); ++i) {
        Eigen::MatrixXd query = x.row(i);

        resultSet.reset(i);
        kdTree.index_->findNeighbors(resultSet, query.data());

        for (int j = 0; j < resultSet.size(); ++j) {
            output(i, j) = resultSet.index(j) + 1;
        }
    }

    return output;
}
