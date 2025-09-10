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

// [[Rcpp::export]]
Rcpp::List vecchia_U_parts_grouped(
    const Eigen::ArrayXi& blockIndices,
    const Eigen::ArrayXi& blockLastIndex,
    const Eigen::ArrayXi& blockLocalResponseIndex,
    const Eigen::ArrayXi& blockLastResponseOfBlock,
    const Eigen::MatrixXd& xWarped,
    const std::string& covarianceFunction,
    const Eigen::VectorXd& sigmaDeviation,
    const Eigen::VectorXd& sigmaSquaredNugget,
    bool parallel
) {
    int nBlocks = blockLastIndex.size();

    Eigen::ArrayXi blockSize(nBlocks);
    Eigen::ArrayXi blockNResponses(nBlocks);
    for (int i = 0; i < nBlocks; i++) {
        if (i == 0) {
            blockSize[i] = blockLastIndex[i];
            blockNResponses[i] = blockLastResponseOfBlock[i];
        } else {
            blockSize[i] = blockLastIndex[i] - blockLastIndex[i - 1];
            blockNResponses[i] = blockLastResponseOfBlock[i] - blockLastResponseOfBlock[i - 1];
        }
    }
    
    Eigen::ArrayXi blockOutputLastIndex(nBlocks);
    for (int i = 0; i < nBlocks; i++) {
        if (i == 0) {
            blockOutputLastIndex[i] = blockSize[i] * blockNResponses[i];
        } else {
            blockOutputLastIndex[i] = blockOutputLastIndex[i - 1] + blockSize[i] * blockNResponses[i];
        }
    }
    Eigen::VectorXi is(blockOutputLastIndex[nBlocks - 1]);
    Eigen::VectorXi js(blockOutputLastIndex[nBlocks - 1]);
    Eigen::VectorXd xs(blockOutputLastIndex[nBlocks - 1]);

    auto fn = [&](int start, int end) {
        for (int k = start; k < end; ++k) {
            Eigen::MatrixXd blockSigma(blockSize[k], blockSize[k]);
            Eigen::LLT<Eigen::MatrixXd> lltSigma(blockSize[k]);
            Eigen::MatrixXd B(blockSize[k], blockNResponses[k]);

            int blockStartIndex = 0;
            int blockResponseStartIndex = 0;
            if (k > 0) {
                blockStartIndex = blockLastIndex[k - 1];
                blockResponseStartIndex = blockLastResponseOfBlock[k - 1];
            }

            B.setZero();
            for (int i = 0; i < blockNResponses[k]; ++i) {
                B(blockLocalResponseIndex[blockResponseStartIndex + i] - 1, i) = 1;
            }


            for (int i = 0; i < blockSize[k]; ++i) {
                int indexI = blockIndices[blockStartIndex + i] - 1;
                for (int j = i + 1; j < blockSize[k]; ++j) {
                    int indexJ = blockIndices[blockStartIndex + j] - 1;
                    double distance = (
                        xWarped.row(indexI) - xWarped.row(indexJ)
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
                        sigmaDeviation[indexI]
                        * sigmaDeviation[indexJ]
                        * correlation
                    );
                    blockSigma(j, i) = blockSigma(i, j);
                }

                blockSigma(i, i) = sigmaDeviation[indexI] * sigmaDeviation[indexI] + sigmaSquaredNugget[indexI];
            }

            lltSigma.compute(blockSigma);

            int blockOutputStartIndex = 0;
            if (k > 0) {
                blockOutputStartIndex = blockOutputLastIndex[k - 1];
            }
            xs.segment(
                blockOutputStartIndex,
                blockSize[k] * blockNResponses[k]
            ) = lltSigma.matrixU().solve(B).reshaped();

            for (int i = 0; i < blockNResponses[k]; ++i) {
                int indexI = blockIndices[blockStartIndex + blockLocalResponseIndex[blockResponseStartIndex + i] - 1];
                for (int j = 0; j < blockSize[k]; ++j) {
                    int indexJ = blockIndices[blockStartIndex + j];

                    is[blockOutputStartIndex + i * blockSize[k] + j] = indexJ;
                    js[blockOutputStartIndex + i * blockSize[k] + j] = indexI;
                }
            }
        }
    };


    #if RCPP_PARALLEL_USE_TBB
    if (parallel) {
        tbb::parallel_for(
            tbb::blocked_range<int>(0, nBlocks),
            [&](const tbb::blocked_range<int>& r) {
                fn(r.begin(), r.end());
            }
        );
    } else {
        fn(0, nBlocks);
    }
#else
    fn(0, nBlocks);
#endif

    Rcpp::List output;
    output["i"] = Rcpp::wrap(is);
    output["j"] = Rcpp::wrap(js);
    output["x"] = Rcpp::wrap(xs);
    return output;
}

