#include <Rcpp.h>

// [[Rcpp::export(name = ".quad_form_diag", rng = false)]]
Rcpp::NumericMatrix quadFormDiag(Rcpp::NumericMatrix A, Rcpp::NumericVector x) {
    if (A.nrow() != A.ncol()) {
        Rcpp::stop("Matrix A must be symmetric");
    }
    int n = A.nrow();
    Rcpp::NumericMatrix B(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B(i, j) = A(i, j) * x[i] * x[j];
        }
    }
    return B;
}

// [[Rcpp::export(name = ".unique_rows", rng = false)]]
Rcpp::NumericMatrix uniqueRows(Rcpp::NumericMatrix input) {
    std::set<std::vector<double>> unique;
    std::vector<double> tmp(input.ncol());
    for (int i = 0; i < input.nrow(); i++) {
        for (int j = 0; j < input.ncol(); j++) {
            tmp[j] = input(i, j);
        }
        unique.insert(tmp);
    }
    Rcpp::NumericMatrix output(unique.size(), input.ncol());
    for (int i = 0; i < output.nrow(); i++) {
        for (int j = 0; j < output.ncol(); j++) {
            output(i, j) = (*std::next(unique.begin(), i))[j];
        }
    }
    Rcpp::IntegerVector indices(input.nrow());
    for (int i = 0; i < input.nrow(); i++) {
        for (int j = 0; j < output.ncol(); j++) {
            tmp[j] = input(i, j);
        }
        indices[i] = std::distance(
            unique.begin(),
            unique.find(tmp)
        ) + 1;
    }
    output.attr("indices") = indices;
    return output;
}
