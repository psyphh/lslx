#include <RcppEigen.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cfloat>
#include <string.h>
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::depends(RcppEigen)]]
class lslxOptimizer {
public:
  std::string algorithm;
  int iter_in_max, iter_out_max, iter_other_max, iter_armijo_max;
  double tol_in, tol_out, tol_other;
  double step_size, armijo;
  double ridge_cov, ridge_hessian;
  bool positive_diag;
  bool response, regularizer;
  
  std::string regularizer_type;
  double lambda, delta;
  int iter_out;
  
  int n_observation;
  Rcpp::List  sample_proportion, saturated_cov, saturated_mean;
  Rcpp::List  saturated_moment_acov;
  
  int n_response, n_factor, n_eta, n_moment, n_group, n_theta;
  
  Rcpp::CharacterVector theta_name;
  Rcpp::LogicalVector theta_is_free, theta_is_pen, theta_is_diag;
  Rcpp::IntegerVector theta_matrice_idx, theta_group_idx;
  Rcpp::IntegerVector theta_left_idx, theta_right_idx, theta_flat_idx;
  
  double baseline_loss_value;
  int baseline_degree_of_freedom;
  
  Eigen::MatrixXd identity_y, identity_eta, identity_theta;  
  Eigen::SparseMatrix<double> identity_y2, duplication_y;
  Eigen::SparseMatrix<double> elimination_y, duplication_eta, commutation_y;
  
  Rcpp::NumericVector theta_start, theta_value, theta_direction;
  Rcpp::IntegerVector theta_est_idx;
  
  Rcpp::List alpha, beta, beta_pinv, psi;
  Rcpp::List mu, sigma, sigma_inv;
  Rcpp::List alpha_derivative, beta_derivative, psi_derivative;
  
  Rcpp::List residual_weight;
  Rcpp::List moment_jacobian;
  
  double loss_value;
  Eigen::MatrixXd loss_gradient;
  Eigen::MatrixXd loss_gradient_diff;
  Eigen::MatrixXd loss_expected_hessian;
  Eigen::MatrixXd loss_observed_hessian;
  Eigen::MatrixXd loss_bfgs_hessian;
  Eigen::MatrixXd loss_bfgs_hessian_inv;
  
  double regularizer_value;
  Eigen::MatrixXd regularizer_gradient;
  
  double objective_value;
  Eigen::MatrixXd objective_gradient;
  
  double objective_gradient_abs_max, objective_hessian_convexity;
  int n_iter_out, n_nonzero_coefficient, degree_of_freedom;
  double robust_degree_of_freedom, scaling_factor;
  
  double aic, aic3, caic;
  double bic, abic, hbic;
  double raic, raic3, rcaic;
  double rbic, rabic, rhbic;
  double rmsea, srmr, cfi, nnfi;
  
  lslxOptimizer(Rcpp::List reduced_data,
                Rcpp::List reduced_model,
                Rcpp::List control,
                Rcpp::List supplied_result);
  
  void set_regularizer(Rcpp::CharacterVector regularizer_type_, double lambda_, double delta_);
  void set_theta_value(Rcpp::NumericVector theta_value_);
  void update_coefficient_matrice();
  void update_implied_moment();
  void update_residual_weight();
  void update_moment_jacobian();
  void update_loss_value();
  void update_loss_gradient();
  void update_loss_gradient_direct();
  void update_loss_expected_hessian();
  void update_loss_observed_hessian();
  void update_loss_bfgs_hessian();
  void update_regularizer_value();
  void update_regularizer_gradient();
  void update_objective_value();
  void update_objective_gradient();
  void update_theta_direction();
  void update_theta_value();
  void update_theta_start();
  void update_numerical_condition();
  void update_information_criterion();
  void update_fit_indice();
  void update_coefficient();
  Rcpp::NumericVector extract_numerical_condition();
  Rcpp::NumericVector extract_information_criterion();
  Rcpp::NumericVector extract_fit_indice();
  Rcpp::NumericVector extract_coefficient();
  Eigen::MatrixXd slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx);
  Eigen::MatrixXd slice_both(Eigen::MatrixXd x, 
                             Rcpp::IntegerVector row_idx, 
                             Rcpp::IntegerVector col_idx);
  Eigen::MatrixXd expand_both(Eigen::MatrixXd x,
                              Rcpp::IntegerVector row_idx, 
                              Rcpp::IntegerVector col_idx,
                              int n_row,
                              int n_col);
  Eigen::MatrixXd vech(Eigen::MatrixXd x);
  Eigen::SparseMatrix<double> create_commutation(int n);
  Eigen::SparseMatrix<double> create_duplication(int n);
  int sign(double x);
};

lslxOptimizer::lslxOptimizer(Rcpp::List reduced_data,
                             Rcpp::List reduced_model,
                             Rcpp::List control,
                             Rcpp::List supplied_result) {
  algorithm = Rcpp::as<std::string>(control["algorithm"]);
  iter_in_max = Rcpp::as<int>(control["iter_in_max"]);
  iter_out_max = Rcpp::as<int>(control["iter_out_max"]);
  iter_other_max =  Rcpp::as<int>(control["iter_other_max"]);
  iter_armijo_max = Rcpp::as<int>(control["iter_armijo_max"]);
  
  tol_in = Rcpp::as<double>(control["tol_in"]);
  tol_out = Rcpp::as<double>(control["tol_out"]);
  tol_other = Rcpp::as<double>(control["tol_other"]);
  
  step_size = Rcpp::as<double>(control["step_size"]);
  ridge_cov = Rcpp::as<double>(control["ridge_cov"]);
  ridge_hessian = Rcpp::as<double>(control["ridge_hessian"]);
  armijo = Rcpp::as<double>(control["armijo"]);
  positive_diag = Rcpp::as<bool>(control["positive_diag"]);
  response = Rcpp::as<bool>(control["response"]);
  regularizer = Rcpp::as<bool>(control["regularizer"]);
  iter_out = -1;
  
  n_response = Rcpp::as<int>(reduced_model["n_response"]);
  n_factor = Rcpp::as<int>(reduced_model["n_factor"]);
  n_eta = Rcpp::as<int>(reduced_model["n_eta"]);
  n_group = Rcpp::as<int>(reduced_model["n_group"]);
  n_moment = Rcpp::as<int>(reduced_model["n_moment"]);
  n_theta = Rcpp::as<int>(reduced_model["n_theta"]);
  
  theta_name = Rcpp::as<CharacterVector>(reduced_model["theta_name"]);
  theta_is_free = Rcpp::as<LogicalVector>(reduced_model["theta_is_free"]);
  theta_is_pen = Rcpp::as<LogicalVector>(reduced_model["theta_is_pen"]);
  theta_is_diag = Rcpp::as<LogicalVector>(reduced_model["theta_is_diag"]);
  theta_matrice_idx = Rcpp::as<IntegerVector>(reduced_model["theta_matrice_idx"]);
  theta_group_idx = Rcpp::as<IntegerVector>(reduced_model["theta_group_idx"]);
  
  theta_left_idx = Rcpp::as<IntegerVector>(reduced_model["theta_left_idx"]) - 1;
  theta_right_idx = Rcpp::as<IntegerVector>(reduced_model["theta_right_idx"]) - 1;
  theta_flat_idx = Rcpp::as<IntegerVector>(reduced_model["theta_flat_idx"]) - 1;
  
  identity_y.resize(n_response, n_response);
  identity_y.setIdentity();
  identity_eta.resize(n_eta, n_eta);
  identity_eta.setIdentity();
  identity_theta.resize(n_theta, n_theta);
  identity_theta.setIdentity();
  identity_y2.resize(n_response * n_response, n_response * n_response);
  identity_y2.setIdentity();
  
  duplication_y  = create_duplication(n_response);
  elimination_y  = duplication_y.transpose() * duplication_y;
  Eigen::SparseLU<SparseMatrix<double> > solver;
  Eigen::SparseMatrix<double> identity_n_cov(n_moment - n_response, 
                                             n_moment - n_response);
  identity_n_cov.setIdentity();
  solver.compute(elimination_y);
  elimination_y = solver.solve(identity_n_cov) * duplication_y.transpose();  
  duplication_eta  = create_duplication(n_eta);
  commutation_y  = create_commutation(n_response);
  
  n_observation = Rcpp::as<int>(reduced_data["n_observation"]);
  sample_proportion = Rcpp::as<List>(reduced_data["sample_proportion"]);
  saturated_cov = Rcpp::as<List>(reduced_data["saturated_cov"]);
  saturated_mean = Rcpp::as<List>(reduced_data["saturated_mean"]);
  saturated_moment_acov = Rcpp::as<List>(reduced_data["saturated_moment_acov"]);
  
  int i;
  for (i = 0; i < n_group; i ++) {
    alpha.push_back(Eigen::MatrixXd::Zero(n_eta, 1));
    beta.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    beta_pinv.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    psi.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    alpha_derivative.push_back(Eigen::MatrixXd::Zero(n_eta, 1));
    beta_derivative.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    psi_derivative.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    mu.push_back(Eigen::MatrixXd::Zero(n_response, 1));
    sigma.push_back(Eigen::MatrixXd::Zero(n_response, n_response));
    sigma_inv.push_back(Eigen::MatrixXd::Zero(n_response, n_response));
    residual_weight.push_back(Eigen::MatrixXd::Zero(n_moment, n_moment));
    moment_jacobian.push_back(Eigen::MatrixXd::Zero(n_moment, n_theta));
  }
  
  baseline_loss_value = Rcpp::as<double>(Rcpp::as<Rcpp::NumericVector>(supplied_result["baseline_model"])["loss_value"]);
  baseline_degree_of_freedom = Rcpp::as<double>(Rcpp::as<Rcpp::NumericVector>(supplied_result["baseline_model"])["degree_of_freedom"]);
  
  theta_start = Rcpp::clone(Rcpp::as<NumericVector>(supplied_result["fitted_start"]));
  theta_value = Rcpp::clone(Rcpp::as<NumericVector>(supplied_result["fitted_start"]));
  theta_direction = Rcpp::rep(0.0, n_theta);
  theta_value.attr("names") = theta_name;
  Rcpp::LogicalVector theta_est_idc = theta_is_pen | theta_is_free;
  for (i = 0; i < n_theta; i++) {
    if (theta_est_idc[i]) {
      theta_est_idx.push_back(i);
    }
  }
  loss_gradient = Eigen::MatrixXd::Zero(n_theta, 1);
  loss_gradient_diff = Eigen::MatrixXd::Zero(n_theta, 1);
  loss_expected_hessian = Eigen::MatrixXd::Identity(n_theta, n_theta);
  loss_observed_hessian = Eigen::MatrixXd::Identity(n_theta, n_theta);
  loss_bfgs_hessian = Eigen::MatrixXd::Identity(n_theta, n_theta);
  loss_bfgs_hessian_inv = Eigen::MatrixXd::Identity(n_theta, n_theta);
  regularizer_gradient = Eigen::MatrixXd::Zero(n_theta, 1);
  objective_gradient= Eigen::MatrixXd::Zero(n_theta, 1);
}

void lslxOptimizer::set_regularizer(Rcpp::CharacterVector regularizer_type_,
                                    double lambda_, 
                                    double delta_) {
  regularizer_type = Rcpp::as<std::string>(regularizer_type_[0]);
  lambda = lambda_;
  delta = delta_;
}

void lslxOptimizer::set_theta_value(Rcpp::NumericVector theta_value_) {
  theta_value = Rcpp::clone(theta_value_);
}

void lslxOptimizer::update_coefficient_matrice() {
  Rcpp::IntegerVector theta_group_idx_unique = Rcpp::sort_unique(theta_group_idx);
  
  Rcpp::IntegerVector theta_left_idx_i0;
  Rcpp::IntegerVector theta_right_idx_i0;
  Rcpp::NumericVector theta_value_i0;
  
  Rcpp::IntegerVector theta_left_idx_ij;
  Rcpp::IntegerVector theta_right_idx_ij;
  Rcpp::NumericVector theta_value_ij;
  
  int theta_left_idx_i0k, theta_right_idx_i0k;
  double theta_value_i0k;
  
  int theta_left_idx_ijk, theta_right_idx_ijk;
  double theta_value_ijk;
  
  int i, j, k;
  for (i = 1; i <= 3; i ++) {
    if (Rcpp::is_true(Rcpp::any(0 == theta_group_idx_unique))) {
      theta_left_idx_i0 =
        theta_left_idx[(theta_matrice_idx == i) &
        (theta_group_idx == 0)];
      theta_right_idx_i0 =
        theta_right_idx[(theta_matrice_idx == i) &
        (theta_group_idx == 0)];
      theta_value_i0 =
        theta_value[(theta_matrice_idx == i) &
        (theta_group_idx == 0)];
    }
    
    for (j = 1; j <= n_group; j ++) {
      if (Rcpp::is_true(Rcpp::any(j == theta_group_idx_unique))) {
        theta_left_idx_ij =
          theta_left_idx[(theta_matrice_idx == i) &
          (theta_group_idx == j)];
        theta_right_idx_ij =
          theta_right_idx[(theta_matrice_idx == i) &
          (theta_group_idx == j)];
        theta_value_ij =
          theta_value[(theta_matrice_idx == i) &
          (theta_group_idx == j)];
        
        switch(i) {
        case 1: {
            Eigen::Map<MatrixXd> alpha_i(Rcpp::as< Eigen::Map <MatrixXd> >(alpha[j - 1]));
            for (k = 0; k < theta_value_ij.size(); k ++) {
              theta_left_idx_ijk = theta_left_idx_ij[k];
              theta_right_idx_ijk = theta_right_idx_ij[k];
              theta_value_ijk = theta_value_ij[k];
              alpha_i(theta_left_idx_ijk, theta_right_idx_ijk) = theta_value_ijk;
            }
            break;
          }
        case 2: {
          Eigen::Map<MatrixXd> beta_i(Rcpp::as< Eigen::Map<MatrixXd> >(beta[j - 1]));
          int k;
          for (k = 0; k < theta_value_ij.size(); k ++) {
            theta_left_idx_ijk = theta_left_idx_ij[k];
            theta_right_idx_ijk = theta_right_idx_ij[k];
            theta_value_ijk = theta_value_ij[k];
            beta_i(theta_left_idx_ijk, theta_right_idx_ijk) = theta_value_ijk;
          }
          break;
        }
        case 3: {
          Eigen::Map<MatrixXd> psi_i(Rcpp::as< Eigen::Map <MatrixXd> >(psi[j - 1]));
          int k;
          for (k = 0; k < theta_value_ij.size(); k ++) {
            theta_left_idx_ijk = theta_left_idx_ij[k];
            theta_right_idx_ijk = theta_right_idx_ij[k];
            theta_value_ijk = theta_value_ij[k];
            psi_i(theta_left_idx_ijk, theta_right_idx_ijk) = theta_value_ijk;
            psi_i(theta_right_idx_ijk, theta_left_idx_ijk) = theta_value_ijk;
          }
          break;
        }
        }
      }
      
      if (Rcpp::is_true(Rcpp::any(0 == theta_group_idx_unique))) {
        switch(i) {
        case 1: {
        Eigen::Map<MatrixXd> alpha_i(Rcpp::as< Eigen::Map <MatrixXd> >(alpha[j - 1]));
        for (k = 0; k < theta_value_i0.size(); k ++) {
          theta_left_idx_i0k = theta_left_idx_i0[k];
          theta_right_idx_i0k = theta_right_idx_i0[k];
          theta_value_i0k = theta_value_i0[k];
          if (Rcpp::is_true(Rcpp::any(j == theta_group_idx_unique))) {
            alpha_i(theta_left_idx_i0k, theta_right_idx_i0k) =
              alpha_i(theta_left_idx_i0k, theta_right_idx_i0k) + theta_value_i0k;
          } else {
            alpha_i(theta_left_idx_i0k, theta_right_idx_i0k) = theta_value_i0k;
          }
        }
        break;
      }
        case 2: {
          Eigen::Map<MatrixXd> beta_i(Rcpp::as< Eigen::Map <MatrixXd> >(beta[j - 1]));
          for (k = 0; k < theta_value_i0.size(); k ++) {
            theta_left_idx_i0k = theta_left_idx_i0[k];
            theta_right_idx_i0k = theta_right_idx_i0[k];
            theta_value_i0k = theta_value_i0[k];
            if (Rcpp::is_true(Rcpp::any(j == theta_group_idx_unique))) {
              beta_i(theta_left_idx_i0k, theta_right_idx_i0k) =
                beta_i(theta_left_idx_i0k, theta_right_idx_i0k) + theta_value_i0k;
            } else {
              beta_i(theta_left_idx_i0k, theta_right_idx_i0k) = theta_value_i0k;
            }
          }
          break;
        }
        case 3: {
          Eigen::Map<MatrixXd> psi_i(Rcpp::as< Eigen::Map <MatrixXd> >(psi[j - 1]));
          for (k = 0; k < theta_value_i0.size(); k ++) {
            theta_left_idx_i0k = theta_left_idx_i0[k];
            theta_right_idx_i0k = theta_right_idx_i0[k];
            theta_value_i0k = theta_value_i0[k];
            if (Rcpp::is_true(Rcpp::any(j == theta_group_idx_unique))) {
              psi_i(theta_left_idx_i0k, theta_right_idx_i0k) =
                psi_i(theta_left_idx_i0k, theta_right_idx_i0k) + theta_value_i0k;
              psi_i(theta_right_idx_i0k, theta_left_idx_i0k) =
                psi_i(theta_left_idx_i0k, theta_right_idx_i0k);
            } else {
              psi_i(theta_left_idx_i0k, theta_right_idx_i0k) = theta_value_i0k;
              psi_i(theta_right_idx_i0k, theta_left_idx_i0k) = theta_value_i0k;
            }
          }
          break;
        }
        }
      }
    }
  }
}



void lslxOptimizer::update_implied_moment() {
  int i;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> alpha_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(alpha[i]));
    Eigen::Map<Eigen::MatrixXd> beta_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta[i]));
    Eigen::Map<Eigen::MatrixXd> psi_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(psi[i]));
    
    Eigen::Map<MatrixXd> beta_pinv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta_pinv[i]));
    Eigen::Map<MatrixXd> mu_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(mu[i]));  
    Eigen::Map<MatrixXd> sigma_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma[i]));
    Eigen::Map<MatrixXd> sigma_inv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma_inv[i]));
    
    beta_pinv_i = (identity_eta - beta_i).inverse();
    mu_i =  beta_pinv_i.topRows(n_response) * alpha_i;
    sigma_i =  beta_pinv_i.topRows(n_response) * psi_i * beta_pinv_i.topRows(n_response).transpose();
    sigma_inv_i = sigma_i.inverse();
  }
}


void lslxOptimizer::update_residual_weight() {
  int i;
  double sample_proportion_i;
  for (i = 0; i < n_group; i++) {
    sample_proportion_i = Rcpp::as<double>(sample_proportion[i]);
    Eigen::Map<MatrixXd> sigma_inv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma_inv[i]));
    Eigen::Map<MatrixXd> residual_weight_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(residual_weight[i]));
    residual_weight_i.block(0, 0, n_response, n_response) = 
      2 * sample_proportion_i * sigma_inv_i;
    residual_weight_i.block(n_response, n_response, 
                            (n_moment - n_response), 
                            (n_moment - n_response)) = 
                              sample_proportion_i * duplication_y.transpose() * 
                              Eigen::kroneckerProduct(sigma_inv_i, sigma_inv_i) * duplication_y;
  }
}


void lslxOptimizer::update_moment_jacobian() {
  Rcpp::IntegerVector theta_group_idx_unique = Rcpp::sort_unique(theta_group_idx);
  Rcpp::IntegerVector theta_flat_idx_j;
  Rcpp::IntegerVector theta_matrice_idx_j;
  Rcpp::IntegerVector theta_flat_idx_jk;
  int n_theta_sum, n_theta_j, n_theta_jk;
  
  int i, j, k;
  for (i = 0; i < n_group; i++) {
    n_theta_sum = 0;
    Eigen::Map<Eigen::MatrixXd> alpha_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(alpha[i]));
    Eigen::Map<Eigen::MatrixXd> beta_pinv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta_pinv[i]));
    Eigen::Map<Eigen::MatrixXd> psi_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(psi[i]));
    Eigen::Map<Eigen::MatrixXd> moment_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(moment_jacobian[i]));
    
    for (j = 0; j < theta_group_idx_unique.size(); j++) {
      n_theta_j = 0;
      theta_flat_idx_j = theta_flat_idx[(theta_group_idx == theta_group_idx_unique[j])];
      theta_matrice_idx_j = theta_matrice_idx[(theta_group_idx == theta_group_idx_unique[j])];
      
      for (k = 1; k <= 3; k++) {
        theta_flat_idx_jk = theta_flat_idx_j[theta_matrice_idx_j == k];
        n_theta_jk = theta_flat_idx_jk.size();
        
        if ((theta_group_idx_unique[j] == 0) | ((theta_group_idx_unique[j] - 1) == i)) {
          if (n_theta_jk > 0) {
            switch(k) {
            case 1: {
          moment_jacobian_i.block(
            0, n_theta_sum + n_theta_j,
            n_response, n_theta_jk) = 
              slice_col(beta_pinv_i.topRows(n_response), theta_flat_idx_jk);
          break;
        }
            case 2: {
              moment_jacobian_i.block(
                0, n_theta_sum + n_theta_j,
                n_response, n_theta_jk) = 
                  slice_col(
                    kroneckerProduct(
                      (alpha_i.transpose() * beta_pinv_i.transpose()),
                      beta_pinv_i.topRows(n_response)),
                      theta_flat_idx_jk);
              
              moment_jacobian_i.block(
                n_response, n_theta_sum + n_theta_j,
                (n_moment - n_response), n_theta_jk) = 
                  slice_col(
                    ((elimination_y * (commutation_y + identity_y2)) * 
                      kroneckerProduct(
                        (beta_pinv_i.topRows(n_response) * psi_i * beta_pinv_i.transpose()), 
                        beta_pinv_i.topRows(n_response))),
                        theta_flat_idx_jk);
              break;
            }
            case 3: {
              moment_jacobian_i.block(
                n_response, n_theta_sum + n_theta_j,
                (n_moment - n_response), n_theta_jk) = 
                  slice_col(
                    (elimination_y * 
                      kroneckerProduct(
                        beta_pinv_i.topRows(n_response), 
                        beta_pinv_i.topRows(n_response)) * duplication_eta),
                        theta_flat_idx_jk);
              break;
            }
            }
          }
        }
        n_theta_j = n_theta_j + n_theta_jk;
      } 
      n_theta_sum = n_theta_sum + n_theta_j;
    }
  }
}

void lslxOptimizer::update_loss_value() {
  loss_value = 0;
  int i;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> saturated_mean_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_mean[i]));
    Eigen::Map<Eigen::MatrixXd> saturated_cov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_cov[i]));
    Eigen::Map<Eigen::MatrixXd> mu_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(mu[i]));
    Eigen::Map<Eigen::MatrixXd> sigma_inv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma_inv[i]));    
    double sample_proportion_i = Rcpp::as<double>(sample_proportion[i]);
    double loss_value_i =
      (saturated_cov_i * sigma_inv_i).diagonal().sum() - 
      log((saturated_cov_i * sigma_inv_i).determinant()) - n_response + 
      ((saturated_mean_i - mu_i).transpose() * sigma_inv_i * (saturated_mean_i - mu_i)).value();
    loss_value += sample_proportion_i * loss_value_i;
  }
}

void lslxOptimizer::update_loss_gradient_direct() {
  Eigen::MatrixXd weight_mu_i, weight_sigma_i;
  double sample_proportion_i;
  loss_gradient_diff = loss_gradient;
  loss_gradient = Eigen::MatrixXd::Zero(n_theta, 1);
  int i, j;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> saturated_mean_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_mean[i]));
    Eigen::Map<Eigen::MatrixXd> saturated_cov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_cov[i]));
    Eigen::Map<Eigen::MatrixXd> mu_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(mu[i]));  
    Eigen::Map<Eigen::MatrixXd> sigma_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma[i]));
    Eigen::Map<Eigen::MatrixXd> sigma_inv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma_inv[i]));
    Eigen::Map<Eigen::MatrixXd> alpha_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(alpha[i]));
    Eigen::Map<Eigen::MatrixXd> beta_pinv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta_pinv[i]));
    Eigen::Map<Eigen::MatrixXd> psi_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(psi[i]));
    Eigen::Map<Eigen::MatrixXd> alpha_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(alpha_derivative[i]));
    Eigen::Map<Eigen::MatrixXd> beta_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta_derivative[i]));
    Eigen::Map<Eigen::MatrixXd> psi_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(psi_derivative[i]));
    weight_mu_i = 2 * (saturated_mean_i - mu_i).transpose() * sigma_inv_i;
    weight_sigma_i = sigma_inv_i * (saturated_cov_i - sigma_i + (saturated_mean_i - mu_i) * (saturated_mean_i - mu_i).transpose()) * sigma_inv_i;
    alpha_derivative_i = - (weight_mu_i * beta_pinv_i.topRows(n_response)).transpose();
    beta_derivative_i = - (
      (beta_pinv_i.topRows(n_response).transpose() * weight_mu_i.transpose() * 
        alpha_i.transpose() * beta_pinv_i.transpose()) +  
        2 * (beta_pinv_i.topRows(n_response).transpose() * weight_sigma_i * 
        beta_pinv_i.topRows(n_response) * psi_i * beta_pinv_i.transpose()));
    psi_derivative_i = - 2 * (beta_pinv_i.topRows(n_response).transpose() * weight_sigma_i * beta_pinv_i.topRows(n_response));
    for (j = 0; j < n_eta; j++) {
      psi_derivative_i(j, j) = 0.5 * psi_derivative_i(j, j);
    } 
  }
  for (i = 0; i < n_theta; i++) {
    if (theta_group_idx[i] == 0) {
      for (j = 0; j < n_group; j++) {
        sample_proportion_i = Rcpp::as<double>(sample_proportion[j]);
        if (theta_matrice_idx[i] == 1) {
          Eigen::Map<Eigen::MatrixXd> alpha_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(alpha_derivative[j]));
          loss_gradient(i, 0) += sample_proportion_i * alpha_derivative_i(theta_left_idx[i], theta_right_idx[i]);
        } else if (theta_matrice_idx[i] == 2) {
          Eigen::Map<Eigen::MatrixXd> beta_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta_derivative[j]));
          loss_gradient(i, 0) += sample_proportion_i * beta_derivative_i(theta_left_idx[i], theta_right_idx[i]);
        } else if (theta_matrice_idx[i] == 3) {
          Eigen::Map<Eigen::MatrixXd> psi_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(psi_derivative[j]));
          loss_gradient(i, 0) += sample_proportion_i * psi_derivative_i(theta_left_idx[i], theta_right_idx[i]);
        } else {}
      }
    } else {
      sample_proportion_i = Rcpp::as<double>(sample_proportion[theta_group_idx[i] - 1]);
      if (theta_matrice_idx[i] == 1) {
        Eigen::Map<Eigen::MatrixXd> alpha_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(alpha_derivative[theta_group_idx[i] - 1]));
        loss_gradient(i, 0) += sample_proportion_i * alpha_derivative_i(theta_left_idx[i], theta_right_idx[i]);
      } else if (theta_matrice_idx[i] == 2) {
        Eigen::Map<Eigen::MatrixXd> beta_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta_derivative[theta_group_idx[i] - 1]));
        loss_gradient(i, 0) += sample_proportion_i * beta_derivative_i(theta_left_idx[i], theta_right_idx[i]);
      } else if (theta_matrice_idx[i] == 3) {
        Eigen::Map<Eigen::MatrixXd> psi_derivative_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(psi_derivative[theta_group_idx[i] - 1]));
        loss_gradient(i, 0) += sample_proportion_i * psi_derivative_i(theta_left_idx[i], theta_right_idx[i]);
      } else {}
    }
  }
  loss_gradient_diff = loss_gradient - loss_gradient_diff;
}



void lslxOptimizer::update_loss_gradient() {
  loss_gradient_diff = loss_gradient;
  loss_gradient = Eigen::MatrixXd::Zero(n_theta, 1);
  Eigen::MatrixXd moment_residual_i = Eigen::MatrixXd::Zero(n_moment, 1);
  int i;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> saturated_mean_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_mean[i]));
    Eigen::Map<Eigen::MatrixXd> saturated_cov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_cov[i]));
    Eigen::Map<Eigen::MatrixXd> mu_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(mu[i]));
    Eigen::Map<Eigen::MatrixXd> sigma_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma[i]));
    Eigen::Map<Eigen::MatrixXd> residual_weight_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(residual_weight[i]));
    Eigen::Map<Eigen::MatrixXd> moment_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(moment_jacobian[i]));
    moment_residual_i.block(
      0, 0, n_response, 1) = (saturated_mean_i - mu_i);
    moment_residual_i.block(
      n_response, 0, (n_moment - n_response), 1) = 
        vech(saturated_cov_i + saturated_mean_i * saturated_mean_i.transpose() - 
        mu_i * saturated_mean_i.transpose() -
        saturated_mean_i * mu_i.transpose() + 
        mu_i * mu_i.transpose() - sigma_i);
    loss_gradient += - moment_jacobian_i.transpose() * residual_weight_i * moment_residual_i;
  }
  loss_gradient_diff = loss_gradient - loss_gradient_diff;
}


void lslxOptimizer::update_loss_expected_hessian() {
  loss_expected_hessian = Eigen::MatrixXd::Zero(n_theta, n_theta);
  int i;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> residual_weight_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(residual_weight[i]));
    Eigen::Map<Eigen::MatrixXd> moment_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(moment_jacobian[i]));
    loss_expected_hessian += moment_jacobian_i.transpose() * residual_weight_i * moment_jacobian_i;
  }
}


void lslxOptimizer::update_loss_observed_hessian() {
  loss_observed_hessian = Eigen::MatrixXd::Zero(n_theta, n_theta);
  Eigen::MatrixXd loss_gradient_0 = loss_gradient;
  int i;
  for (i = 0; i < n_theta; i++) {
    theta_value = Rcpp::clone(theta_start);
    theta_value[i] = theta_value[i] + tol_other;
    update_coefficient_matrice();
    update_implied_moment();
    update_loss_gradient_direct(); 
    loss_observed_hessian.col(i) = (loss_gradient - loss_gradient_0) / tol_other;
  }
  loss_observed_hessian = (loss_observed_hessian + loss_observed_hessian.transpose()) / 2.0;
  theta_value = Rcpp::clone(theta_start);
  loss_gradient = loss_gradient_0;
}


void lslxOptimizer::update_loss_bfgs_hessian() {
  Eigen::VectorXd theta_diff = 
    Rcpp::as<Eigen::VectorXd>(theta_value) - Rcpp::as<Eigen::VectorXd>(theta_start);
  double rho = 1.0 / (loss_gradient_diff.transpose() * theta_diff).value();
  int i;
  if (iter_out <= 0) {
    loss_bfgs_hessian = loss_expected_hessian;
    loss_bfgs_hessian_inv = 
      expand_both(slice_both(loss_expected_hessian, theta_est_idx, theta_est_idx).inverse(),
                  theta_est_idx, theta_est_idx,
                  n_theta, n_theta);
  } else {
    for (i = 0; i < n_theta; i++) {
      if (!(theta_is_free[i] || theta_is_pen[i])) {
        loss_gradient_diff(i, 0) = 0;
      }
    }
    loss_bfgs_hessian = loss_bfgs_hessian - 
      (loss_bfgs_hessian * theta_diff * theta_diff.transpose() * loss_bfgs_hessian) /
        (theta_diff.transpose() * loss_bfgs_hessian * theta_diff).value() +
          rho * (loss_gradient_diff * loss_gradient_diff.transpose());
    loss_bfgs_hessian_inv = (identity_theta - rho * theta_diff * loss_gradient_diff.transpose()) * 
      loss_bfgs_hessian_inv * (identity_theta - rho * loss_gradient_diff * theta_diff.transpose()) + 
      rho * theta_diff * theta_diff.transpose();
  }
}

void lslxOptimizer::update_regularizer_value() {
  regularizer_value = 0.0;
  int i;
  double regularizer_value_i;
  if (lambda > DBL_EPSILON) {
    for (i = 0; i < n_theta; i++) {
      if (theta_is_pen[i]) {
        if (std::fabs(theta_value[i]) < (lambda * delta)) {
          if (std::fabs(theta_value[i]) < DBL_EPSILON) {
            regularizer_value_i = 0.0;
          } else {
            regularizer_value_i = 
              lambda * (std::fabs(theta_value[i]) - std::pow(theta_value[i], 2) / (2.0 * lambda * delta));
          }
        } else {
          regularizer_value_i = (std::pow(lambda, 2) * delta) / 2.0;
        }
      } else {
        regularizer_value_i = 0;
      }
      regularizer_value += regularizer_value_i;
    }
  } else {
  }
}

void lslxOptimizer::update_regularizer_gradient() {
  int i;
  if (lambda > DBL_EPSILON) {
    for (i = 0; i < n_theta; i++) {
      if (theta_is_pen[i]) {
        if ((theta_value[i] <= (lambda * delta)) & (theta_value[i] > DBL_EPSILON)) {
          regularizer_gradient(i, 0) = lambda - (theta_value[i] / delta);
        } else if ((- theta_value[i] <= (lambda * delta)) & (theta_value[i] < - DBL_EPSILON)) {
          regularizer_gradient(i, 0) = - lambda - (theta_value[i] / delta);
        } else if ((theta_value[i] > (lambda * delta)) | ((- theta_value[i]) > (lambda * delta))) {
          regularizer_gradient(i, 0) = 0;
        } else {
          regularizer_gradient(i, 0) = lambda;
        }
      } else {
        regularizer_gradient(i, 0) = 0;
      }
    }
  } else {
    regularizer_gradient = Eigen::MatrixXd::Zero(n_theta, 1);
  }
}

void lslxOptimizer::update_objective_value() {
  objective_value = loss_value + regularizer_value;
}

void lslxOptimizer::update_objective_gradient() {
  int i;
  for (i = 0; i < n_theta; i++) {
    if (std::fabs(theta_value[i]) > DBL_EPSILON) {
      objective_gradient(i, 0) = loss_gradient(i, 0) + regularizer_gradient(i, 0);
    } else {
      objective_gradient(i, 0) = sign(loss_gradient(i, 0)) * 
        std::max((std::fabs(loss_gradient(i, 0)) - lambda), 0.0);
    }
  }
}

void lslxOptimizer::update_theta_direction() {
  theta_direction = Rcpp::rep(0.0, n_theta);
  Rcpp::NumericVector z = Rcpp::rep(0.0, n_theta);
  Eigen::MatrixXd g, h;
  double z_r, z_l;
  int i, j;
  if (regularizer) {
    g = loss_gradient;
    if (algorithm == "bfgs") {
      h = loss_bfgs_hessian;
    } else if (algorithm == "fisher") {
      h = loss_expected_hessian;
    } else {}
    for (i = 0; i < n_theta; i++) {
      h(i, i) = h(i, i) + ridge_hessian;
    }
    for (i = 0; i < iter_in_max; i++) {
      for (j = 0; j < n_theta; j++) {
        Eigen::Map<Eigen::VectorXd> d(Rcpp::as< Eigen::Map <Eigen::VectorXd> >(theta_direction));
        double g_ij = g(j, 0) + (h * d)(j, 0);
        double h_ij = h(j, j);
        if (lambda > DBL_EPSILON) {
          if (theta_is_free[j]) {
            z[j] = (-g_ij / h_ij);
          } else if (theta_is_pen[j]) {
            z_r = ((theta_value[j] + theta_direction[j]) / delta - g_ij - lambda) / (h_ij - (1.0 / delta));
            z_l = ((theta_value[j] + theta_direction[j]) / delta - g_ij + lambda) / (h_ij - (1.0 / delta));
            if (z_r >= - (theta_value[j] + theta_direction[j])) {
              if (z_r >= (lambda * delta - (theta_value[j] + theta_direction[j]))) {
                z[j] = (-g_ij / h_ij);
              } else {
                z[j] = z_r;
              }
            } else if (z_l <= -(theta_value[j] + theta_direction[j])) {
              if (z_l <= (-lambda * delta - (theta_value[j] + theta_direction[j]))) {
                z[j] = (-g_ij / h_ij);
              } else {
                z[j] = z_l;
              }
            } else {
              z[j] = - (theta_value[j] + theta_direction[j]);
            }
          } else {
            z[j] = 0;
          }
        } else {
          if (theta_is_free[j] | theta_is_pen[j]) {
            z[j] = (-g_ij / h_ij);
          } else {
            z[j] = 0;
          }
        }
        theta_direction[j] = theta_direction[j] + z[j];
      }
      if ((h.diagonal().array() * Rcpp::as<Eigen::VectorXd>(z).array().abs()).maxCoeff() < tol_in) {
        break;
      }
    }
  } else {
    g = loss_gradient;
    if (algorithm == "bfgs") {
      h = loss_bfgs_hessian_inv;
    } else if (algorithm == "fisher") {
      h = expand_both(slice_both(
        loss_expected_hessian, theta_est_idx, theta_est_idx).inverse(),
        theta_est_idx, theta_est_idx,
        n_theta, n_theta);
    } else {}
    theta_direction = - h * g;
  }
  double theta_direction_norm = Rcpp::as<Eigen::VectorXd>(theta_direction).norm();
  if (theta_direction_norm > 1) {
    theta_direction = theta_direction / theta_direction_norm;
  }
}


void lslxOptimizer::update_theta_value() {
  Rcpp::IntegerVector theta_group_idx_unique = Rcpp::sort_unique(theta_group_idx);
  double objective_value_old = objective_value;
  double regularizer_value_old = regularizer_value;
  double step_size_i;
  double regularizer_value_0 = regularizer_value;
  int i;
  for (i = 0; i < iter_armijo_max; i++) {
    step_size_i = std::pow(step_size, i);
    theta_value = theta_start + step_size_i * theta_direction;
    if (positive_diag) {
      if (!Rcpp::is_true(Rcpp::any(0 == theta_group_idx_unique))) {
        theta_value = 
          Rcpp::ifelse((theta_value < 0) & (theta_is_diag), ridge_cov, theta_value);
      } else {
        theta_value = 
          Rcpp::ifelse((theta_value < 0) & theta_is_diag & (theta_group_idx == 0), ridge_cov, theta_value);
      }
    }
    update_coefficient_matrice();
    update_implied_moment();
    update_loss_value();
    update_regularizer_value();
    update_objective_value();
    if (i == 0) {
      regularizer_value_0 = regularizer_value;
    }
    if ((objective_value - objective_value_old) <= 
        ((armijo * step_size_i) * 
        ((Rcpp::as<Eigen::VectorXd>(theta_direction).transpose() * loss_gradient).value() +
        (regularizer_value_0 - regularizer_value_old))) ) {
      break;
    }
  }
}


void lslxOptimizer::update_theta_start() {
  theta_start = Rcpp::clone(theta_value);
}


void lslxOptimizer::update_coefficient() {
  Rcpp::NumericVector objective_gradient_abs(n_theta);
  if (iter_out == -1) {
    update_coefficient_matrice();
    update_implied_moment();
    update_loss_value();
    update_residual_weight();
    update_moment_jacobian();
    update_loss_gradient();
    update_loss_expected_hessian();
  } 
  if (algorithm == "bfgs") {
    update_loss_bfgs_hessian();
  } 
  update_regularizer_value();
  update_objective_value();
  
  int i;  
  for (iter_out = 1; iter_out <= iter_out_max; iter_out++) {
    update_theta_direction();
    update_theta_value();
    if (algorithm == "bfgs") {
      update_loss_gradient_direct();
      update_loss_bfgs_hessian();
    } else if (algorithm == "fisher") {
      update_residual_weight();
      update_moment_jacobian();
      update_loss_gradient();
      update_loss_expected_hessian();
    } else {}
    update_regularizer_gradient();
    update_objective_gradient();
    update_theta_start();
    for (i = 0; i < n_theta; i++) {
      if (theta_is_free[i] | theta_is_pen[i]) {
        objective_gradient_abs[i] = std::fabs(objective_gradient(i, 0));
      } else {
        objective_gradient_abs[i] = - INFINITY;
      }
    }
    objective_gradient_abs_max = Rcpp::max(objective_gradient_abs);
    n_iter_out = iter_out;
    if ((objective_gradient_abs_max < tol_out) | (iter_out == iter_out_max)) {
      iter_out = 0;
      break;
    }
  }
}

void lslxOptimizer::update_numerical_condition() {
  Rcpp::NumericVector objective_hessian_diagonal(n_theta);
  Rcpp::IntegerVector idx_is_effective(0);
  Eigen::MatrixXd loss_hessian;
  Eigen::MatrixXd saturated_moment_acov_matrix, residual_weight_matrix, moment_jacobian_matrix;
  if (algorithm == "bfgs") {
    loss_hessian = loss_bfgs_hessian;
  } else if (algorithm == "fisher") {
    loss_hessian = loss_expected_hessian;    
  } else{}
  int i;
  for (i = 0; i < n_theta; i++) {
    if (theta_is_free[i]) {
      objective_hessian_diagonal[i] = loss_hessian(i, i) + ridge_hessian;
      idx_is_effective.push_back(i);
    } else if (theta_is_pen[i]) {
      objective_hessian_diagonal[i] = loss_hessian(i, i) + ridge_hessian - (1 / delta);
      if (std::abs(theta_value[i]) > DBL_EPSILON) {
        idx_is_effective.push_back(i);
      } else {
      }
    } else {
      objective_hessian_diagonal[i] = INFINITY;
    }
  }
  objective_hessian_convexity = Rcpp::min(objective_hessian_diagonal);
  n_nonzero_coefficient = idx_is_effective.size(); 
  degree_of_freedom = n_group * n_moment - n_nonzero_coefficient;
  if (response) {
    update_moment_jacobian();
    update_residual_weight();
    if (idx_is_effective.size() > 0) {
      saturated_moment_acov_matrix = Eigen::MatrixXd::Zero(n_group * n_moment, n_group * n_moment);
      residual_weight_matrix = Eigen::MatrixXd::Zero(n_group * n_moment, n_group * n_moment); 
      moment_jacobian_matrix = Eigen::MatrixXd::Zero(n_group * n_moment, n_theta);
      for (i = 0; i < n_group; i++) {
        Eigen::Map<Eigen::MatrixXd> saturated_moment_acov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_moment_acov[i]));
        Eigen::Map<Eigen::MatrixXd> residual_weight_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(residual_weight[i]));
        Eigen::Map<Eigen::MatrixXd> moment_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(moment_jacobian[i]));
        saturated_moment_acov_matrix.block(
          i * n_moment, i * n_moment,
          n_moment, n_moment) = saturated_moment_acov_i;
        residual_weight_matrix.block(
          i * n_moment, i * n_moment,
          n_moment, n_moment) = residual_weight_i;
        moment_jacobian_matrix.block(
          i * n_moment, 0,
          n_moment, n_theta) = moment_jacobian_i;
      }
      moment_jacobian_matrix = slice_col(moment_jacobian_matrix, idx_is_effective);
      robust_degree_of_freedom = 0.5 * double(n_observation) * (saturated_moment_acov_matrix * 
        (residual_weight_matrix - (residual_weight_matrix * moment_jacobian_matrix) *
        (moment_jacobian_matrix.transpose() * residual_weight_matrix * moment_jacobian_matrix).inverse() *
        (moment_jacobian_matrix.transpose() * residual_weight_matrix))).diagonal().sum();
      if (degree_of_freedom > 0) {
        scaling_factor = robust_degree_of_freedom / degree_of_freedom;
      } else {
        scaling_factor = NAN;
      }
    } else {
      scaling_factor = NAN;      
    }
  } else {
    robust_degree_of_freedom = NAN;
    scaling_factor = NAN;
  }
}




void lslxOptimizer::update_information_criterion() {
  aic = loss_value - (2.0 / double(n_observation)) * double(degree_of_freedom);
  aic3 = loss_value - (3.0 / double(n_observation)) * double(degree_of_freedom);
  caic = loss_value - ((1 + std::log(double(n_observation))) / double(n_observation)) * double(degree_of_freedom);
  
  bic = loss_value - (std::log(double(n_observation)) / double(n_observation)) * double(degree_of_freedom);
  abic = loss_value - (std::log((double(n_observation) + 2.0) / 24.0) / double(n_observation)) * double(degree_of_freedom);
  hbic = loss_value - (std::log(double(n_observation) / (2.0 * 3.1415926)) / double(n_observation)) * double(degree_of_freedom);
  
  raic = loss_value - (2.0 / double(n_observation)) * double(robust_degree_of_freedom);
  raic3 = loss_value - (3.0 / double(n_observation)) * double(robust_degree_of_freedom);
  rcaic = loss_value - ((1 + std::log(double(n_observation))) / double(n_observation)) * double(robust_degree_of_freedom);
  
  rbic = loss_value - (std::log(double(n_observation)) / double(n_observation)) * double(robust_degree_of_freedom);
  rabic = loss_value - (std::log((double(n_observation) + 2.0) / 24.0) / double(n_observation)) * double(robust_degree_of_freedom);
  rhbic = loss_value - (std::log(double(n_observation) / (2.0 * 3.1415926)) / double(n_observation)) * double(robust_degree_of_freedom);
  
}


void lslxOptimizer::update_fit_indice() {
  if ((degree_of_freedom == 0) & (loss_value > std::sqrt(DBL_EPSILON))) {
    rmsea = NAN;
  } else {
    if (loss_value < std::sqrt(DBL_EPSILON)) {
      rmsea = 0;
    } else {
      rmsea = std::sqrt(n_group * std::max(((loss_value / double(degree_of_freedom)) - 
        (1 / double(n_observation))), 0.0)); 
    }
  }
  
  double cfi_num = std::max((double(n_observation) * loss_value - double(degree_of_freedom)), 0.0);
  double cfi_den = std::max(std::max(double(n_observation) * loss_value - double(degree_of_freedom),
                                     double(n_observation) * baseline_loss_value - double(baseline_degree_of_freedom)), 0.0);
  if ((cfi_num < std::sqrt(DBL_EPSILON)) & (cfi_den < std::sqrt(DBL_EPSILON))) {
    cfi = NAN;
  } else {
    if (cfi_num < DBL_EPSILON) {
      cfi = 1;
    } else {
      cfi = 1 - cfi_num / cfi_den;
    }
  }
  double nnfi_0 = (double(n_observation) * baseline_loss_value) / double(baseline_degree_of_freedom);
  double nnfi_1;
  if ((loss_value > std::sqrt(DBL_EPSILON)) & (degree_of_freedom == 0)) {
    nnfi = NAN;
  } else {
    if (loss_value < std::sqrt(DBL_EPSILON)) {
      nnfi_1 = 0;
    } else {
      nnfi_1 = (double(n_observation) * loss_value) / double(degree_of_freedom);
    }
    nnfi = (nnfi_0 - nnfi_1) / (nnfi_0 - 1.0);
    nnfi = std::min(nnfi, 1.0);
  }
  
  srmr = 0;
  int i, j, k;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> saturated_mean_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_mean[i]));
    Eigen::Map<Eigen::MatrixXd> saturated_cov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_cov[i]));
    Eigen::Map<Eigen::MatrixXd> mu_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(mu[i]));
    Eigen::Map<Eigen::MatrixXd> sigma_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma[i]));    
    double sample_proportion_i = Rcpp::as<double>(sample_proportion[i]);
    double residual_sigma_i = 0;
    double residual_mu_i = 0;
    for (j = 0; j < n_response; j++) {
      for (k = j; k < n_response; k++) {
        residual_sigma_i += std::pow((saturated_cov_i(j, k) - sigma_i(j, k)), 2.0) / (sigma_i.diagonal()(j) * sigma_i.diagonal()(k));
      }
      residual_mu_i += std::pow((saturated_mean_i(j) - mu_i(j)), 2.0) / (sigma_i.diagonal()(j));
    }
    double srmr_sigma_i =  residual_sigma_i / (double(n_response) * (double(n_response) + 1.0) / 2.0); 
    double srmr_mu_i = residual_mu_i / double(n_response); 
    double srmr_i = sample_proportion_i * std::sqrt(srmr_sigma_i + srmr_mu_i);
    srmr += srmr_i;
  }
}

Rcpp::NumericVector lslxOptimizer::extract_numerical_condition() {
  Rcpp::NumericVector numerical_condition = 
    Rcpp::NumericVector::create(
      _["lambda"] = lambda,
      _["delta"] = delta,
      _["objective_value"] = objective_value,
      _["objective_gradient_abs_max"] = objective_gradient_abs_max,
      _["objective_hessian_convexity"] = objective_hessian_convexity,
      _["n_iter_out"] = n_iter_out,
      _["loss_value"] = loss_value,
      _["n_nonzero_coefficient"] = n_nonzero_coefficient,
      _["degree_of_freedom"] = degree_of_freedom,
      _["robust_degree_of_freedom"] = robust_degree_of_freedom,
      _["scaling_factor"] = scaling_factor);
  return Rcpp::clone(numerical_condition);
}

Rcpp::NumericVector lslxOptimizer::extract_information_criterion() {
  Rcpp::NumericVector information_criterion = 
    Rcpp::NumericVector::create(
      _["aic"] = aic,
      _["aic3"] = aic3,
      _["caic"] = caic,
      _["bic"] = bic,
      _["abic"] = abic,
      _["hbic"] = hbic,
      _["raic"] = raic,
      _["raic3"] = raic3,
      _["rcaic"] = rcaic,
      _["rbic"] = rbic,
      _["rabic"] = rabic,
      _["rhbic"] = rhbic);
  return Rcpp::clone(information_criterion);
}

Rcpp::NumericVector lslxOptimizer::extract_fit_indice() {
  Rcpp::NumericVector fit_indice = 
    Rcpp::NumericVector::create(
      _["rmsea"] = rmsea,
      _["cfi"] = cfi,
      _["nnfi"] = nnfi,
      _["srmr"] = srmr);
  return Rcpp::clone(fit_indice);
}

Rcpp::NumericVector lslxOptimizer::extract_coefficient() {
  return Rcpp::clone(theta_value);
}



Eigen::MatrixXd lslxOptimizer::slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx) {
  Eigen::MatrixXd y(x.rows(), col_idx.size());
  int i;
  for (i = 0; i < col_idx.size(); i++) {
    y.col(i) = x.col(col_idx[i]);
  }
  return(y);
}



Eigen::MatrixXd lslxOptimizer::slice_both(Eigen::MatrixXd x, 
                                          Rcpp::IntegerVector row_idx, 
                                          Rcpp::IntegerVector col_idx) {
  Eigen::MatrixXd y(row_idx.size(), col_idx.size());
  int i, j;
  for (i = 0; i < row_idx.size(); i++) {
    for (j = 0; j < col_idx.size(); j++) {
      y(i, j) = x(row_idx[i], col_idx[j]);
    }
  }
  return y;
}

Eigen::MatrixXd lslxOptimizer::expand_both(Eigen::MatrixXd x, 
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
  return y;
}



Eigen::MatrixXd lslxOptimizer::vech(Eigen::MatrixXd x) {
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
  return y;
}

Eigen::SparseMatrix<double> lslxOptimizer::create_commutation(int n) {
  int n2 = n * n;
  Eigen::SparseMatrix<double> commutation(n2, n2);
  int i, j, row_idx, col_idx;
  for (i = 0; i < n2; i++) {
    row_idx = i % n;
    col_idx = i / n;
    j = n * row_idx + col_idx;
    commutation.insert(i, j) = 1;
  }
  return commutation;
}

Eigen::SparseMatrix<double> lslxOptimizer::create_duplication(int n) {
  SparseMatrix<double> duplication(n * n, (n * (n + 1)) / 2);
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
        duplication.insert(idx, col_idx) = 1;
      }
      idx += 1;
    }
  }
  return duplication;
}





int lslxOptimizer::sign(double x) {
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


Eigen::MatrixXd slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx) {
  Eigen::MatrixXd y(x.rows(), col_idx.size());
  int i;
  for (i = 0; i < col_idx.size(); i++) {
    y.col(i) = x.col(col_idx[i]);
  }
  return(y);
}


Eigen::MatrixXd slice_row(Eigen::MatrixXd x, Rcpp::IntegerVector row_idx) {
  Eigen::MatrixXd y(row_idx.size(), x.cols());
  int i;
  for (i = 0; i < row_idx.size(); i++) {
    y.row(i) = x.row(row_idx[i]);
  }
  return(y);
}


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


Eigen::MatrixXd expand_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx, int n_col) {
  Eigen::MatrixXd y;
  y = Eigen::MatrixXd::Zero(x.rows(), n_col);
  int i;
  for (i = 0; i < col_idx.size(); i++) {
    y.col(col_idx[i]) = x.col(i);
  }
  return(y);
}



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


Eigen::SparseMatrix<double> create_duplication(int n) {
  SparseMatrix<double> duplication(n * n, (n * (n + 1)) / 2);
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
        duplication.insert(idx, col_idx) = 1;
      }
      idx += 1;
    }
  }
  return duplication;
}



// [[Rcpp::export]]
void compute_regularized_path_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result,
    Rcpp::List fitted_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Rcpp::NumericVector lambda_grid = Rcpp::as<Rcpp::NumericVector>(control["lambda_grid"]);
  Rcpp::NumericVector delta_grid = Rcpp::as<Rcpp::NumericVector>(control["delta_grid"]);
  Rcpp::List numerical_condition = Rcpp::as<Rcpp::List>(fitted_result["numerical_condition"]);
  Rcpp::List information_criterion = Rcpp::as<Rcpp::List>(fitted_result["information_criterion"]);
  Rcpp::List fit_indice = Rcpp::as<Rcpp::List>(fitted_result["fit_indice"]);
  Rcpp::List coefficient = Rcpp::as<Rcpp::List>(fitted_result["coefficient"]);
  
  int i, j, idx;
  for (i = 0; i < lambda_grid.size(); i++) {
    for (j = 0; j < delta_grid.size(); j++) {
      optimizer.set_regularizer(
        Rcpp::as< Rcpp::CharacterVector >(control["penalty_method"]), 
        lambda_grid[i], delta_grid[j]);
      optimizer.update_coefficient();
      optimizer.update_numerical_condition();
      optimizer.update_information_criterion();
      optimizer.update_fit_indice();
      idx = i * delta_grid.size() + j;
      coefficient[idx] = optimizer.extract_coefficient();
      numerical_condition[idx] = optimizer.extract_numerical_condition();
      information_criterion[idx] = optimizer.extract_information_criterion();
      fit_indice[idx] = optimizer.extract_fit_indice();
    }
  }
}

// [[Rcpp::export]]
Rcpp::List compute_coefficient_matrice_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Rcpp::List coefficient_matrice;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  
  coefficient_matrice = 
    Rcpp::List::create(Rcpp::Named("alpha") = optimizer.alpha,
                       Rcpp::Named("beta") = optimizer.beta,
                       Rcpp::Named("psi") = optimizer.psi);
  
  return Rcpp::wrap(coefficient_matrice);
}

// [[Rcpp::export]]
Rcpp::List compute_implied_cov_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Rcpp::List implied_cov;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  implied_cov = optimizer.sigma;
  return Rcpp::wrap(implied_cov);
}

// [[Rcpp::export]]
Rcpp::List compute_implied_mean_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Rcpp::List implied_mean;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  implied_mean = optimizer.mu;
  return Rcpp::wrap(implied_mean);
}



// [[Rcpp::export]]
Rcpp::NumericMatrix compute_moment_jacobian_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Eigen::MatrixXd moment_jacobian = 
    Eigen::MatrixXd::Zero(optimizer.n_group * optimizer.n_moment, 
                          optimizer.n_theta);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  optimizer.update_moment_jacobian();
  int i;
  for (i = 0; i < optimizer.n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> moment_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.moment_jacobian[i]));
    moment_jacobian.block(
      i * optimizer.n_moment, 0,
      optimizer.n_moment, optimizer.n_theta) = moment_jacobian_i;
  }
  return Rcpp::wrap(moment_jacobian);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_expected_fisher_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd expected_fisher;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  optimizer.update_residual_weight();
  optimizer.update_moment_jacobian();
  optimizer.update_loss_expected_hessian();
  expected_fisher = 0.5 * optimizer.loss_expected_hessian;
  return Rcpp::wrap(expected_fisher);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_observed_fisher_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd observed_fisher;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_theta_start();
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  optimizer.update_residual_weight();
  optimizer.update_moment_jacobian();
  optimizer.update_loss_gradient();
  optimizer.update_loss_observed_hessian();
  observed_fisher = 0.5 * optimizer.loss_observed_hessian;
  return Rcpp::wrap(observed_fisher);
}



// [[Rcpp::export]]
Rcpp::NumericMatrix compute_score_acov_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Eigen::MatrixXd score_acov = 
    Eigen::MatrixXd::Zero(optimizer.n_theta, 
                          optimizer.n_theta);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  optimizer.update_residual_weight();
  optimizer.update_moment_jacobian();
  
  int i;
  for (i = 0; i < optimizer.n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> moment_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.moment_jacobian[i]));
    Eigen::Map<Eigen::MatrixXd> residual_weight_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.residual_weight[i]));
    Eigen::Map<Eigen::MatrixXd> saturated_moment_acov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.saturated_moment_acov[i]));
    score_acov += 0.25 *
      moment_jacobian_i.transpose() * residual_weight_i * 
      saturated_moment_acov_i * residual_weight_i * moment_jacobian_i;
  }
  return Rcpp::wrap(score_acov);
}



// [[Rcpp::export]]
double compute_loss_value_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  double loss_value;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  optimizer.update_loss_value();
  loss_value = optimizer.loss_value;
  return loss_value;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix compute_loss_gradient_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd loss_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  optimizer.update_residual_weight();
  optimizer.update_moment_jacobian();
  optimizer.update_loss_gradient();
  loss_gradient = optimizer.loss_gradient;
  return Rcpp::wrap(loss_gradient);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_loss_gradient_direct_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd loss_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  optimizer.update_loss_gradient_direct();
  loss_gradient = optimizer.loss_gradient;
  return Rcpp::wrap(loss_gradient);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_regularizer_gradient_cpp(
    Rcpp::NumericVector theta_value,
    double lambda,
    double delta,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd regularizer_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.set_regularizer(Rcpp::as<Rcpp::CharacterVector>(control["penalty_method"]), lambda, delta);
  optimizer.update_regularizer_gradient();
  regularizer_gradient = optimizer.regularizer_gradient;
  return Rcpp::wrap(regularizer_gradient);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_objective_gradient_cpp(
    Rcpp::NumericVector theta_value,
    double lambda,
    double delta,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd objective_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.set_regularizer(Rcpp::as<Rcpp::CharacterVector>(control["penalty_method"]), lambda, delta);
  
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  optimizer.update_loss_gradient_direct();
  optimizer.update_regularizer_gradient();
  optimizer.update_objective_gradient();
  objective_gradient = optimizer.objective_gradient;
  return Rcpp::wrap(objective_gradient);
}


// [[Rcpp::export]]
void compute_saturated_moment_cpp(
    Rcpp::List y_obs,
    Rcpp::List w,
    Rcpp::List m_idx,
    Rcpp::List saturated_mean,
    Rcpp::List saturated_cov,
    int iter_other_max,
    double tol_other) {
  Rcpp::List y_obs_i;
  Rcpp::List w_i;
  Rcpp::List m_idx_i;
  Eigen::MatrixXd saturated_mean_i, saturated_cov_i;
  Rcpp::IntegerVector m_idx_ik;
  Eigen::MatrixXd saturated_cov_ik, saturated_cov_ik_inv;
  Eigen::MatrixXd e_sum_i, c_sum_i;
  Eigen::MatrixXd a_ik, b_ik;
  Eigen::MatrixXd y_ik, y_ik_w;
  double moment_diff_max;
  
  int i, j, k;
  int n_group = y_obs.size();
  for (i = 0; i < n_group; i++) {
    y_obs_i = Rcpp::as<List>(y_obs[i]);
    m_idx_i = Rcpp::as<List>(m_idx[i]);
    w_i = Rcpp::as<List>(w[i]);
    for (j = 0; j < iter_other_max; j++) {
      saturated_cov_i = Rcpp::as< Eigen::MatrixXd >(saturated_cov[i]);
      saturated_mean_i = Rcpp::as< Eigen::MatrixXd >(saturated_mean[i]);
      e_sum_i = Eigen::MatrixXd::Zero(saturated_mean_i.rows(), 1);
      c_sum_i = Eigen::MatrixXd::Zero(saturated_cov_i.rows(), saturated_cov_i.cols());
      for (k = 0; k < y_obs_i.size(); k++) {
        Eigen::Map<Eigen::MatrixXd> y_obs_ik(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(y_obs_i[k]));
        Eigen::Map<Eigen::VectorXd> w_ik(Rcpp::as< Eigen::Map <Eigen::VectorXd> >(w_i[k]));
        m_idx_ik = Rcpp::as< Rcpp::IntegerVector >(m_idx_i[k]);
        saturated_cov_ik = slice_both(saturated_cov_i, m_idx_ik, m_idx_ik);
        saturated_cov_ik_inv = saturated_cov_ik.inverse();
        
        a_ik = saturated_mean_i.transpose() - 
          slice_row(saturated_mean_i, m_idx_ik).transpose() * 
          saturated_cov_ik_inv * slice_row(saturated_cov_i, m_idx_ik);
        b_ik = saturated_cov_ik_inv * slice_row(saturated_cov_i, m_idx_ik);
        y_ik = y_obs_ik * b_ik + Eigen::MatrixXd::Ones(y_obs_ik.rows(), 1)  * a_ik;
        y_ik_w = (y_ik.array().colwise() * w_ik.array()).matrix();
        e_sum_i += y_ik_w.colwise().sum().transpose();
        c_sum_i += w_ik.sum() * saturated_cov_i -
          w_ik.sum() * slice_col(saturated_cov_i, m_idx_ik) * b_ik +
          y_ik_w.transpose() * y_ik;
      }
      saturated_mean[i] = e_sum_i;
      saturated_cov[i] = c_sum_i - e_sum_i * e_sum_i.transpose();
      moment_diff_max = 
        std::max((Rcpp::as< Eigen::MatrixXd >(saturated_mean[i]) - saturated_mean_i).array().abs().maxCoeff(),
                 (Rcpp::as< Eigen::MatrixXd >(saturated_cov[i]) - saturated_cov_i).array().abs().maxCoeff());
      if (moment_diff_max < tol_other) {
        break;
      }
    }
  }
}




// [[Rcpp::export]]
void compute_saturated_moment_acov_response_cpp(
    Rcpp::List y_obs,
    Rcpp::List w,
    Rcpp::List m_idx,
    Rcpp::List m2_idx,
    Rcpp::List saturated_mean,
    Rcpp::List saturated_cov,
    Rcpp::List saturated_moment_acov) {
  int n_response_i, n_moment_i, sample_size_i;
  Rcpp::List y_obs_i;
  Rcpp::List w_i;
  Rcpp::List m_idx_i;
  Rcpp::List m2_idx_i;
  Eigen::MatrixXd saturated_mean_i, saturated_cov_i;
  Eigen::SparseMatrix<double> duplication_i;
  Eigen::MatrixXd score2_sum_i;
  Eigen::MatrixXd hessian_sum_i, hessian_sum_i_inv;
  Eigen::MatrixXd saturated_mean_ij, saturated_cov_ij, saturated_cov_ij_vech, saturated_cov_ij_inv;
  Eigen::MatrixXd y_obs_ij;
  Eigen::VectorXd w_ij;
  Eigen::MatrixXd yc_obs_ij, yc_obs_ij_w;
  Eigen::MatrixXd yc2c_obs_ij, yc2c_obs_ijk;
  Eigen::MatrixXd score_ij, score_ij_w;
  Eigen::MatrixXd saturated_moment_acov_i;
  Rcpp::IntegerVector m_idx_ij, m2_idx_ij;
  int n_response_ij, n_moment_ij, sample_size_ij;
  Eigen::SparseMatrix<double> duplication_ij;
  
  int i, j, k;
  int n_group = y_obs.size();
  for (i = 0; i < n_group; i++) {
    sample_size_i = 0;
    y_obs_i = Rcpp::as<List>(y_obs[i]);
    m_idx_i = Rcpp::as<List>(m_idx[i]);
    m2_idx_i = Rcpp::as<List>(m2_idx[i]);
    w_i = Rcpp::as<List>(w[i]);
    saturated_cov_i = Rcpp::as< Eigen::MatrixXd >(saturated_cov[i]);
    saturated_mean_i = Rcpp::as< Eigen::MatrixXd >(saturated_mean[i]);
    n_response_i = saturated_mean_i.rows();
    n_moment_i = (n_response_i * (n_response_i + 3)) / 2;
    score2_sum_i = Eigen::MatrixXd::Zero(n_moment_i, n_moment_i);
    hessian_sum_i = Eigen::MatrixXd::Zero(n_moment_i, n_moment_i);
    
    duplication_i = create_duplication(n_response_i);
    if ((y_obs_i.size() == 1) & ((Rcpp::as< Rcpp::IntegerVector >(m_idx_i[0])).size() == n_response_i)) {
      y_obs_ij = Rcpp::as<Eigen::MatrixXd>(y_obs_i[0]);
      w_ij = Rcpp::as<Eigen::VectorXd>(w_i[0]);
      sample_size_ij = y_obs_ij.rows();
      sample_size_i += sample_size_ij;
      saturated_mean_ij = saturated_mean_i;
      saturated_cov_ij = saturated_cov_i;
      saturated_cov_ij_vech = vech(saturated_cov_ij); 
      n_response_ij = saturated_mean_ij.rows();
      n_moment_ij = (n_response_ij * (n_response_ij + 3)) / 2;
      
      yc_obs_ij = y_obs_ij - 
        Eigen::MatrixXd::Ones(sample_size_ij, 1) * saturated_mean_ij.transpose();
      yc2c_obs_ij.resize(sample_size_ij, n_moment_ij - n_response_ij);
      
      for (k = 0; k < sample_size_ij; k++) {
        yc2c_obs_ijk = vech(yc_obs_ij.row(k).transpose() * yc_obs_ij.row(k));
        yc2c_obs_ij.row(k) = (yc2c_obs_ijk - saturated_cov_ij_vech).transpose();
      }
      saturated_cov_ij_inv = saturated_cov_ij.inverse();
      duplication_ij = create_duplication(n_response_ij);
      
      score_ij.resize(sample_size_ij, n_moment_i);
      score_ij.leftCols(n_response_i) = yc_obs_ij;
      score_ij.rightCols((n_moment_i - n_response_i)) = yc2c_obs_ij;
      
      score_ij_w = (score_ij.array().colwise() * w_ij.array()).matrix();
      score2_sum_i += score_ij_w.transpose() * score_ij;
      saturated_moment_acov_i = score2_sum_i / double(sample_size_i);
    } else {
      for (j = 0; j < y_obs_i.size(); j++) {
        y_obs_ij = Rcpp::as<Eigen::MatrixXd>(y_obs_i[j]);
        w_ij = Rcpp::as<Eigen::VectorXd>(w_i[j]);
        m_idx_ij = Rcpp::as< Rcpp::IntegerVector >(m_idx_i[j]);
        m2_idx_ij = Rcpp::as< Rcpp::IntegerVector >(m2_idx_i[j]);
        sample_size_ij = y_obs_ij.rows();
        sample_size_i += sample_size_ij;
        saturated_mean_ij = slice_row(saturated_mean_i, m_idx_ij);
        saturated_cov_ij = slice_both(saturated_cov_i, m_idx_ij, m_idx_ij);
        saturated_cov_ij_vech = vech(saturated_cov_ij); 
        n_response_ij = saturated_mean_ij.rows();
        n_moment_ij = (n_response_ij * (n_response_ij + 3)) / 2;
        
        yc_obs_ij = y_obs_ij - 
          Eigen::MatrixXd::Ones(sample_size_ij, 1) * saturated_mean_ij.transpose();
        yc2c_obs_ij.resize(sample_size_ij, n_moment_ij - n_response_ij);
        
        for (k = 0; k < sample_size_ij; k++) {
          yc2c_obs_ijk = vech(yc_obs_ij.row(k).transpose() * yc_obs_ij.row(k));
          yc2c_obs_ij.row(k) = (yc2c_obs_ijk - saturated_cov_ij_vech).transpose();
        }
        saturated_cov_ij_inv = saturated_cov_ij.inverse();
        duplication_ij = create_duplication(n_response_ij);
        
        score_ij.resize(sample_size_ij, n_moment_i);
        score_ij.leftCols(n_response_i) = 
          expand_col((yc_obs_ij * saturated_cov_ij_inv), m_idx_ij, n_response_i);
        score_ij.rightCols((n_moment_i - n_response_i)) = 
          expand_col((0.5 * yc2c_obs_ij * duplication_ij.transpose() *
          Eigen::kroneckerProduct(saturated_cov_ij_inv, saturated_cov_ij_inv) * duplication_ij), 
          m2_idx_ij, (n_moment_i - n_response_i));
        
        score_ij_w = (score_ij.array().colwise() * w_ij.array()).matrix();
        score2_sum_i += score_ij_w.transpose() * score_ij;
        
        yc_obs_ij = expand_col(yc_obs_ij, m_idx_ij, n_response_i);
        yc_obs_ij_w = (yc_obs_ij.array().colwise() * w_ij.array()).matrix();
        yc2c_obs_ij = expand_col(yc2c_obs_ij, m2_idx_ij, (n_moment_i - n_response_i));
        saturated_cov_ij_inv = expand_both(saturated_cov_ij_inv, m_idx_ij, m_idx_ij,
                                           n_response_i, n_response_i);
        
        hessian_sum_i.block(0, 0, n_response_i, n_response_i) +=
          w_ij.sum() * saturated_cov_ij_inv;
        hessian_sum_i.block(n_response_i, n_response_i, 
                            (n_moment_i - n_response_i), 
                            (n_moment_i - n_response_i)) += 
                              duplication_i.transpose() * 
                              (Eigen::kroneckerProduct(saturated_cov_ij_inv,
                                                       (saturated_cov_ij_inv * 
                                                         (yc_obs_ij_w.transpose() * yc_obs_ij) * 
                                                         saturated_cov_ij_inv - 0.5 * w_ij.sum() * 
                                                         saturated_cov_ij_inv))) * duplication_i;
        hessian_sum_i.block(0, n_response_i, 
                            n_response_i, 
                            (n_moment_i - n_response_i)) += 
                              (Eigen::kroneckerProduct(saturated_cov_ij_inv, 
                                                       yc_obs_ij_w.colwise().sum() * saturated_cov_ij_inv)) * duplication_i;
        hessian_sum_i.block(n_response_i, 0, 
                            (n_moment_i - n_response_i), 
                            n_response_i) = 
                              hessian_sum_i.block(0, n_response_i, 
                                                  n_response_i, 
                                                  (n_moment_i - n_response_i)).transpose(); 
      }
      hessian_sum_i_inv = hessian_sum_i.inverse();
      saturated_moment_acov_i = (hessian_sum_i_inv * score2_sum_i * hessian_sum_i_inv) / double(sample_size_i);
    }
    saturated_moment_acov[i] = saturated_moment_acov_i;
  }
}


// [[Rcpp::export]]
void compute_saturated_moment_acov_moment_cpp(
    int n_observation,
    Rcpp::List sample_proportion,
    Rcpp::List saturated_cov,
    Rcpp::List saturated_moment_acov) {
  Eigen::MatrixXd saturated_cov_i_inv, saturated_moment_acov_i;
  Eigen::SparseMatrix<double> duplication_i;
  double sample_proportion_i;
  int n_response_i, n_moment_i;
  int i;
  for (i = 0; i < saturated_cov.size(); i++) {
    Eigen::Map<MatrixXd> saturated_cov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_cov[i]));
    sample_proportion_i = Rcpp::as<double>(sample_proportion[i]);
    saturated_cov_i_inv = saturated_cov_i.inverse();
    n_response_i = saturated_cov_i.cols();
    n_moment_i = n_response_i * (n_response_i + 3) / 2;
    duplication_i = create_duplication(n_response_i);
    saturated_moment_acov_i = Eigen::MatrixXd::Zero(n_moment_i, n_moment_i);
    saturated_moment_acov_i.block(0, 0, n_response_i, n_response_i) = saturated_cov_i_inv;
    saturated_moment_acov_i.block(n_response_i, n_response_i,
                                  n_moment_i - n_response_i, n_moment_i - n_response_i) =
                                    0.5 * duplication_i.transpose() * 
                                    Eigen::kroneckerProduct(saturated_cov_i_inv, saturated_cov_i_inv) *
                                    duplication_i;
    saturated_moment_acov_i = saturated_moment_acov_i.inverse() / 
      (sample_proportion_i *double(n_observation));
    saturated_moment_acov[i] = saturated_moment_acov_i;
  }
}
