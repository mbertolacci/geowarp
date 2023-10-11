#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

// [[Rcpp::export(name = ".deduplicate_parents", rng = false)]]
Rcpp::IntegerMatrix deduplicateParents(Rcpp::IntegerMatrix parents) {
    using Rcpp::IntegerVector;
    using Rcpp::_;

    Rcpp::IntegerMatrix output(parents.nrow(), parents.ncol());
    for (int i = 0; i < parents.nrow(); i++) {
        IntegerVector parentsI = Rcpp::sort_unique(Rcpp::na_omit(
            parents(i, _)
        ));
        for (int j = 0; j < parents.ncol(); j++) {
            output(i, j) = (
                j < parentsI.size()
                ? parentsI[parentsI.size() - j - 1]
                : NA_INTEGER
            );
        }
    }
    return output;
}

// [[Rcpp::export]]
Rcpp::NumericVector vecchia_U_x_parts(
    const Eigen::MatrixXi& parents,
    const Eigen::MatrixXd& xWarped,
    const std::string& covarianceFunction,
    const Eigen::VectorXd& sigmaDeviation,
    const Eigen::VectorXd& sigmaSquaredNugget,
    bool parallel
) {
    int nRows = parents.rows();
    int blockSize = parents.cols();
    Eigen::VectorXd x(nRows * blockSize);

    auto fn = [&](int start, int end) {
        Eigen::MatrixXd blockSigma(blockSize, blockSize);
        Eigen::LLT<Eigen::MatrixXd> lltSigma(blockSize);
        Eigen::VectorXd b(blockSize);
        b.head(blockSize - 1).setZero();
        b(blockSize - 1) = 1;
        for (int k = start; k < end; ++k) {
            for (int i = 0; i < blockSize; ++i) {
                for (int j = 0; j < blockSize; ++j) {
                    double distance = (
                        xWarped.row(parents(k, i)) - xWarped.row(parents(k, j))
                    ).norm();

                    double correlation;
                    if (covarianceFunction == "exponential") {
                        correlation = exp(-1.0 * distance);
                    } else {
                        correlation = (
                            (1.0 + std::sqrt(3.0) * distance)
                            * std::exp(-std::sqrt(3.0) * distance)
                        );
                    }
                    blockSigma(i, j) = (
                        sigmaDeviation[parents(k, i)]
                        * sigmaDeviation[parents(k, j)]
                        * correlation
                    );
                }

                blockSigma(i, i) += sigmaSquaredNugget[parents(k, i)];
            }

            lltSigma.compute(blockSigma);

            x.segment(
                k * blockSize,
                blockSize
            ) = lltSigma.matrixU().solve(b);
        }
    };

#if RCPP_PARALLEL_USE_TBB
    if (parallel) {
        tbb::parallel_for(
            tbb::blocked_range<int>(0, nRows),
            [&](const tbb::blocked_range<int>& r) {
                fn(r.begin(), r.end());
            }
        );
    } else {
        fn(0, nRows);
    }
#else
    fn(0, nRows);
#endif

    return Rcpp::wrap(x);
}
