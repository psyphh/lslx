// utility functions for package lslx
// written by Po-Hsien Huang psyphh@gmail.com

#include <RcppEigen.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cfloat>
#include <string.h>
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::depends(RcppEigen)]]
// slice columns
Eigen::MatrixXd slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx);

// slice rows
Eigen::MatrixXd slice_row(Eigen::MatrixXd x, Rcpp::IntegerVector row_idx);

// slice both rows and columns
Eigen::MatrixXd slice_both(Eigen::MatrixXd x, 
                           Rcpp::IntegerVector row_idx, 
                           Rcpp::IntegerVector col_idx);

// expand columns
Eigen::MatrixXd expand_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx, int n_col);

// expand both rows and columns
Eigen::MatrixXd expand_both(Eigen::MatrixXd x, 
                            Rcpp::IntegerVector row_idx, 
                            Rcpp::IntegerVector col_idx,
                            int n_row,
                            int n_col);

// vech operator
Eigen::MatrixXd vech(Eigen::MatrixXd x);

// vech operator for only non-diagonal elements
Eigen::MatrixXd vech_small(Eigen::MatrixXd x);

// method for creating commutation matrix
Eigen::SparseMatrix<double> create_commutation(int n);

// create duplication matrix
Eigen::SparseMatrix<double> create_duplication(int n);

// method for which function
Rcpp::IntegerVector which(Rcpp::LogicalVector x);

// method for sign function
int sign(double x);