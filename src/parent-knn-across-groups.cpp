#include <Rcpp.h>
#include <RcppEigen.h>

#include "nanoflann.hpp"

template <
    typename _DistanceType, typename _IndexType = size_t,
    typename _CountType = size_t
>
class GroupParentKNNResultSet {
public:
    using DistanceType = _DistanceType;
    using IndexType    = _IndexType;
    using CountType    = _CountType;

private:
    IndexType target;
    std::vector<IndexType> indices;
    std::vector<DistanceType> dists;
    std::vector<int> groups;
    CountType capacity;
    CountType count;

public:
    explicit GroupParentKNNResultSet(CountType capacity_, const std::vector<int>& groups_)
        : indices(capacity_),
          dists(capacity_),
          groups(groups_),
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

    bool addPoint(DistanceType dist, IndexType index) {
        if (index > target || groups[index] == groups[target]) {
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

// [[Rcpp::export(name = ".parent_knn_across_groups", rng = false)]]
Rcpp::IntegerMatrix parent_knn_across_groups(
    const Eigen::MatrixXd& x,
    Rcpp::IntegerVector groups,
    int n_parents,
    int leaf_size = 40
) {
    typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> KDTreeType;

    if (groups.size() != x.rows()) {
        Rcpp::stop("groups must have the same number of entries as x has rows");
    }

    KDTreeType kdTree(x.cols(), x, leaf_size);
    Rcpp::IntegerMatrix output(x.rows(), n_parents + 1);
    output.fill(Rcpp::IntegerVector::get_na());

    GroupParentKNNResultSet<double> resultSet(
        n_parents,
        Rcpp::as<std::vector<int>>(groups)
    );

    for (int i = 0; i < x.rows(); ++i) {
        Eigen::MatrixXd query = x.row(i);

        resultSet.reset(i);
        kdTree.index_->findNeighbors(resultSet, query.data());

        output(i, 0) = i + 1;
        for (int j = 0; j < resultSet.size(); ++j) {
            output(i, j + 1) = resultSet.index(j) + 1;
        }
    }

    return output;
}
