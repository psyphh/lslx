// utility functions for package lslx
// written by Po-Hsien Huang psyphh@gmail.com

#include "utility-function.h"

// [[Rcpp::depends(RcppEigen)]]
// slice columns
Eigen::MatrixXd slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx) {
  Eigen::MatrixXd y(x.rows(), col_idx.size());
  int i;
  for (i = 0; i < col_idx.size(); i++) {
    y.col(i) = x.col(col_idx[i]);
  }
  return(y);
}

// slice rows
Eigen::MatrixXd slice_row(Eigen::MatrixXd x, Rcpp::IntegerVector row_idx) {
  Eigen::MatrixXd y(row_idx.size(), x.cols());
  int i;
  for (i = 0; i < row_idx.size(); i++) {
    y.row(i) = x.row(row_idx[i]);
  }
  return(y);
}

// slice both rows and columns
Eigen::MatrixXd slice_both(Eigen::MatrixXd x, 
                           Rcpp::IntegerVector row_idx, 
                           Rcpp::IntegerVector col_idx) {
  Eigen::MatrixXd y(row_idx.size(), col_idx.size());
  int i, j;
  for (i = 0; i < row_idx.size(); i++) {
    for (j = 0; j < col_idx.size(); j++) {
      y(i, j) = x(row_idx[i], col_idx[j]);
    }
  }
  return(y);
}

// expand columns
Eigen::MatrixXd expand_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx, int n_col) {
  Eigen::MatrixXd y;
  y = Eigen::MatrixXd::Zero(x.rows(), n_col);
  int i;
  for (i = 0; i < col_idx.size(); i++) {
    y.col(col_idx[i]) = x.col(i);
  }
  return(y);
}

// expand both rows and columns
Eigen::MatrixXd expand_both(Eigen::MatrixXd x, 
                            Rcpp::IntegerVector row_idx, 
                            Rcpp::IntegerVector col_idx,
                            int n_row,
                            int n_col) {
  Eigen::MatrixXd y;
  y = Eigen::MatrixXd::Zero(n_row, n_col);
  int i, j;
  for (i = 0; i < row_idx.size(); i++) {
    for (j = 0; j < col_idx.size(); j++) {
      y(row_idx[i], col_idx[j]) = x(i, j);
    }
  }
  return(y);
}

// vec operator
Eigen::MatrixXd vec(Eigen::MatrixXd x) {
  int n_col = x.cols();
  Eigen::MatrixXd y(n_col * n_col, 1);
  int idx = 0;
  int i, j;
  for (i = 0; i < n_col; i ++ ) {
    for (j = 0; j < n_col; j ++ ) {
      y(idx, 0) = x(j, i);
      idx += 1;
    }
  }
  return(y);
}

// vech operator
Eigen::MatrixXd vech(Eigen::MatrixXd x) {
  int n_col = x.cols();
  Eigen::MatrixXd y((n_col * (n_col + 1)) / 2, 1);
  int idx = 0;
  int i, j;
  for (i = 0; i < n_col; i ++ ) {
    for (j = i; j < n_col; j ++ ) {
      y(idx, 0) = x(j, i);
      idx += 1;
    }
  }
  return(y);
}

// vech operator for only non-diagonal elements
Eigen::MatrixXd vech_small(Eigen::MatrixXd x) {
  int n_col = x.cols();
  Eigen::MatrixXd y((n_col * (n_col - 1)) / 2, 1);
  int idx = 0;
  int i, j;
  for (i = 0; i < (n_col - 1); i ++ ) {
    for (j = (i + 1); j < n_col; j ++ ) {
      y(idx, 0) = x(j, i);
      idx += 1;
    }
  }
  return(y);
}

// method for creating commutation matrix
Eigen::MatrixXd create_commutation(int n) {
  int n2 = n * n;
  Eigen::MatrixXd commutation;
  commutation = Eigen::MatrixXd::Zero(n2, n2);
  int i, j, row_idx, col_idx;
  for (i = 0; i < n2; i++) {
    row_idx = i % n;
    col_idx = i / n;
    j = n * row_idx + col_idx;
    commutation(i, j) = 1.0;
  }
  return commutation;
}


// create duplication matrix
Eigen::MatrixXd create_duplication(int n) {
  Eigen::MatrixXd duplication;
  duplication = Eigen::MatrixXd::Zero(n * n, (n * (n + 1)) / 2);
  int i, j, idx, row_idx, col_idx;
  idx = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j >= i) {
        row_idx = n * i + j - i * (i + 1) / 2;
        col_idx = n * i + j - i * (i + 1) / 2;
      } else {
        row_idx = n * j + i - j * (j + 1) / 2;
        col_idx = n * j + i - j * (j + 1) / 2;
      }
      if (row_idx == col_idx) {
        duplication(idx, col_idx) = 1.0;
      }
      idx += 1;
    }
  }
  return duplication;
}


// method for which function
Rcpp::IntegerVector which(Rcpp::LogicalVector x) {
  Rcpp::IntegerVector y = Rcpp::seq(0, x.size()-1);
  return y[x];
}

// method for sign function
int sign(double x) {
  int y;
  if (x > DBL_EPSILON) {
    y = 1;
  } else if (-x > DBL_EPSILON) {
    y = -1;
  } else {
    y = 0;
  }
  return(y);
}

