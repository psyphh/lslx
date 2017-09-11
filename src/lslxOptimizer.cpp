#include <RcppEigen.h>
#include <limits>
#include <cmath>
#include <cfloat>
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::depends(RcppEigen)]]

class lslxOptimizer {
public:
  char regularizer_type;
  double lambda, gamma;
  int iter_in_max, iter_out_max, iter_armijo_max;
  double tol_in, tol_out;
  double step_size, armijo;
  double ridge_cov, ridge_hessian;
  bool positive_diag;
  
  int total_sample_size;
  Rcpp::List  sample_proportion, saturated_cov, saturated_mean;
  Eigen::MatrixXd  saturated_moment_acov;
  
  int n_response, n_factor, n_eta, n_moment, n_group;
  
  Rcpp::CharacterVector theta_name;
  Rcpp::LogicalVector theta_is_free, theta_is_pen, theta_is_diag;
  Rcpp::IntegerVector theta_matrice_idx, theta_group_idx;
  Rcpp::IntegerVector theta_left_idx, theta_right_idx, theta_flat_idx;
  
  double loss_value_baseline;
  int degree_of_freedom_baseline;

  Eigen::MatrixXd identity_y, identity_eta;  
  Eigen::SparseMatrix<double> identity_y2, duplication_y;
  Eigen::SparseMatrix<double> elimination_y, duplication_eta, commutation_y;
  
  Rcpp::NumericVector theta_start, theta_value, theta_direction;
  
  Rcpp::List alpha, beta, beta_pinv, psi;
  Rcpp::List mu, sigma, sigma_inv;
  
  Eigen::MatrixXd weight_normal;
  Eigen::MatrixXd moment_gradient;
  
  double loss_value;
  Eigen::VectorXd loss_gradient;
  Eigen::MatrixXd loss_expected_hessian;
  
  double regularizer_value;
  Eigen::VectorXd regularizer_gradient;
  
  double objective_value;
  Eigen::VectorXd objective_gradient;
  
  double objective_gradient_abs_max, objective_hessian_convexity;
  int n_iter_out, n_nonzero_coefficient, degree_of_freedom;
  
  double aic, aic3, caic;
  double bic, abic, hbic;
  double rmsea, srmr, cfi, nnfi;
  
  lslxOptimizer(Rcpp::List reduced_data,
                Rcpp::List reduced_model,
                Rcpp::List control);
  
  void set_regularizer(char regularizer_type, double lambda_, double gamma_);
  void set_theta_value(Rcpp::NumericVector theta_value_);
  void update_coefficient_matrice();
  void update_implied_moment();
  void update_weight_normal();
  void update_moment_gradient();
  void update_loss_value();
  void update_loss_gradient();
  void update_loss_expected_hessian();
  void update_regularizer_value();
  void update_regularizer_gradient();
  void update_objective_value();
  void update_objective_gradient();
  void update_theta_direction();
  void update_theta_value();
  void update_theta_start();
  void update_numerical_condition();
  void update_goodness_of_fit();
  void optimize_theta_value();
  Rcpp::NumericVector extract_numerical_condition();
  Rcpp::NumericVector extract_theta_value();
  Rcpp::NumericVector extract_goodness_of_fit();
  Eigen::MatrixXd slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx);
  Eigen::MatrixXd vech(Eigen::MatrixXd x);
  int sign(double x);
};

lslxOptimizer::lslxOptimizer(Rcpp::List reduced_data,
                             Rcpp::List reduced_model,
                             Rcpp::List control) {

  iter_in_max = Rcpp::as<int>(control["iter_in_max"]);
  iter_out_max = Rcpp::as<int>(control["iter_out_max"]);
  iter_armijo_max = Rcpp::as<int>(control["iter_armijo_max"]);
  
  tol_in = Rcpp::as<double>(control["tol_in"]);
  tol_out = Rcpp::as<double>(control["tol_out"]);
  
  step_size = Rcpp::as<double>(control["step_size"]);
  ridge_cov = Rcpp::as<double>(control["ridge_cov"]);
  ridge_hessian = Rcpp::as<double>(control["ridge_hessian"]);
  armijo = Rcpp::as<double>(control["armijo"]);
  positive_diag = Rcpp::as<bool>(control["positive_diag"]);
  
  total_sample_size = Rcpp::as<int>(reduced_data["total_sample_size"]);
  sample_proportion = Rcpp::as<List>(reduced_data["sample_proportion"]);
  saturated_cov = Rcpp::as<List>(reduced_data["saturated_cov"]);
  saturated_mean = Rcpp::as<List>(reduced_data["saturated_mean"]);
  
  n_response = Rcpp::as<int>(reduced_model["n_response"]);
  n_factor = Rcpp::as<int>(reduced_model["n_factor"]);
  n_eta = Rcpp::as<int>(reduced_model["n_eta"]);
  n_group = Rcpp::as<int>(reduced_model["n_group"]);
  n_moment = Rcpp::as<int>(reduced_model["n_moment"]);
  
  theta_name = Rcpp::as<CharacterVector>(reduced_model["theta_name"]);
  theta_is_free = Rcpp::as<LogicalVector>(reduced_model["theta_is_free"]);
  theta_is_pen = Rcpp::as<LogicalVector>(reduced_model["theta_is_pen"]);
  theta_is_diag = Rcpp::as<LogicalVector>(reduced_model["theta_is_diag"]);
  theta_matrice_idx = Rcpp::as<IntegerVector>(reduced_model["theta_matrice_idx"]);
  theta_group_idx = Rcpp::as<IntegerVector>(reduced_model["theta_group_idx"]);
  
  theta_left_idx = Rcpp::as<IntegerVector>(reduced_model["theta_left_idx"]) - 1;
  theta_right_idx = Rcpp::as<IntegerVector>(reduced_model["theta_right_idx"]) - 1;
  theta_flat_idx = Rcpp::as<IntegerVector>(reduced_model["theta_flat_idx"]) - 1;
  
  int i;
  for (i = 0; i < n_group; i ++) {
    alpha.push_back(Eigen::MatrixXd::Zero(n_eta, 1));
    beta.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    beta_pinv.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    psi.push_back(Eigen::MatrixXd::Zero(n_eta, n_eta));
    mu.push_back(Eigen::MatrixXd::Zero(n_response, 1));
    sigma.push_back(Eigen::MatrixXd::Zero(n_response, n_response));
    sigma_inv.push_back(Eigen::MatrixXd::Zero(n_response, n_response));
  }
  
  loss_value_baseline = Rcpp::as<double>(reduced_model["loss_value_baseline"]);
  degree_of_freedom_baseline = Rcpp::as<int>(reduced_model["degree_of_freedom_baseline"]);
  
  identity_y  = Rcpp::as<Eigen::MatrixXd >(reduced_model["identity_y"]);
  identity_eta  = Rcpp::as<Eigen::MatrixXd >(reduced_model["identity_eta"]);
  
  identity_y2  = Rcpp::as<Eigen::SparseMatrix<double> >(reduced_model["identity_y2"]);
  duplication_y  = Rcpp::as<Eigen::SparseMatrix<double> >(reduced_model["duplication_y"]);
  elimination_y  = Rcpp::as<Eigen::SparseMatrix<double> >(reduced_model["elimination_y"]);
  duplication_eta  = Rcpp::as<Eigen::SparseMatrix<double> >(reduced_model["duplication_eta"]);
  commutation_y  = Rcpp::as<Eigen::SparseMatrix<double> >(reduced_model["commutation_y"]);
  
  theta_start = Rcpp::clone(Rcpp::as<NumericVector>(reduced_model["fitted_start"]));
  theta_value = Rcpp::clone(Rcpp::as<NumericVector>(reduced_model["fitted_start"]));
  theta_direction = Rcpp::rep(0.0, theta_name.size());
  theta_value.attr("names") = theta_name;
  
  moment_gradient = Eigen::MatrixXd::Zero(n_group * n_moment, theta_name.size());
  weight_normal = Eigen::MatrixXd::Zero(n_group * n_moment, n_group * n_moment);
  
  loss_gradient = Eigen::MatrixXd::Zero(theta_name.size(), 1);
  loss_expected_hessian = Eigen::MatrixXd::Zero(theta_name.size(), theta_name.size());
  regularizer_gradient = Eigen::MatrixXd::Zero(theta_name.size(), 1);
  objective_gradient= Eigen::MatrixXd::Zero(theta_name.size(), 1);
  
}

void lslxOptimizer::set_regularizer(char regularizer_type_,
                                    double lambda_, 
                                    double gamma_) {
  regularizer_type = regularizer_type_;
  lambda = lambda_;
  gamma = gamma_;
}

void lslxOptimizer::set_theta_value(Rcpp::NumericVector theta_value_) {
  theta_value = theta_value_;
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

void lslxOptimizer::update_weight_normal() {
  int i;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<MatrixXd> sigma_inv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma_inv[i]));
    double sample_proportion_i = Rcpp::as<double>(sample_proportion[i]);
    weight_normal.block(i * n_moment, 
                        i * n_moment, 
                        n_response, n_response) = 
                          2 * sample_proportion_i * sigma_inv_i;
    weight_normal.block((i * n_moment) + n_response, 
                        (i * n_moment) + n_response, 
                        (n_moment - n_response), 
                        (n_moment - n_response)) = 
                          sample_proportion_i * duplication_y.transpose() * 
                          kroneckerProduct(sigma_inv_i, sigma_inv_i) * duplication_y;
  }
}


void lslxOptimizer::update_moment_gradient() {
  Rcpp::IntegerVector theta_group_idx_unique = Rcpp::sort_unique(theta_group_idx);
  
  Rcpp::IntegerVector theta_flat_idx_j;
  Rcpp::IntegerVector theta_matrice_idx_j;
  Rcpp::IntegerVector theta_flat_idx_jk;
  int n_theta_sum, n_theta_j, n_theta_jk;
  
  int i, j, k;
  for (i = 0; i < n_group; i++) {
    n_theta_sum = 0;
    Eigen::Map<Eigen::MatrixXd> alpha_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(alpha[i]));
    Eigen::Map<MatrixXd> beta_pinv_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(beta_pinv[i]));
    Eigen::Map<Eigen::MatrixXd> psi_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(psi[i]));
    
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
          moment_gradient.block(
            (i * n_moment), 
            n_theta_sum + n_theta_j,
            n_response,
            n_theta_jk) = 
              slice_col(beta_pinv_i.topRows(n_response), theta_flat_idx_jk);
          break;
        }
            case 2: {
              moment_gradient.block(
                (i * n_moment), 
                n_theta_sum + n_theta_j,
                n_response,
                n_theta_jk) = 
                  slice_col(
                    kroneckerProduct(
                      (alpha_i.transpose() * beta_pinv_i.transpose()),
                      beta_pinv_i.topRows(n_response)),
                      theta_flat_idx_jk);
              
              moment_gradient.block(
                (i * n_moment) + n_response, 
                n_theta_sum + n_theta_j,
                (n_moment - n_response),
                n_theta_jk) = 
                  slice_col(
                    ((elimination_y * (commutation_y + identity_y2)) * 
                      kroneckerProduct(
                        (beta_pinv_i.topRows(n_response) * psi_i * beta_pinv_i.transpose()), 
                        beta_pinv_i.topRows(n_response))),
                        theta_flat_idx_jk);
              break;
            }
            case 3: {
              moment_gradient.block(
                (i * n_moment) + n_response, 
                n_theta_sum + n_theta_j,
                (n_moment - n_response),
                n_theta_jk) = 
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


void lslxOptimizer::update_loss_gradient() {
  Eigen::MatrixXd moment_residual = Eigen::MatrixXd::Zero(n_group * n_moment, 1);
  int i;
  for (i = 0; i < n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> saturated_mean_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_mean[i]));
    Eigen::Map<Eigen::MatrixXd> saturated_cov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_cov[i]));
    Eigen::Map<Eigen::MatrixXd> mu_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(mu[i]));
    Eigen::Map<Eigen::MatrixXd> sigma_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(sigma[i]));
    moment_residual.block(
      i * n_moment, 0, 
      n_response, 1) = (saturated_mean_i - mu_i);
    moment_residual.block(
      i * n_moment + n_response, 0, 
      (n_moment - n_response), 1) = 
        vech(saturated_cov_i + saturated_mean_i * saturated_mean_i.transpose() - 
        mu_i * saturated_mean_i.transpose() -
        saturated_mean_i * mu_i.transpose() + 
        mu_i * mu_i.transpose() - sigma_i);
  }
  loss_gradient = - moment_gradient.transpose() * weight_normal * moment_residual;
}


void lslxOptimizer::update_loss_expected_hessian() {
  loss_expected_hessian = moment_gradient.transpose() * weight_normal * moment_gradient;
}

void lslxOptimizer::update_regularizer_value() {
  regularizer_value = 0;
  int i;
  double regularizer_value_i;
  if (lambda > DBL_EPSILON) {
    for (i = 0; i < theta_name.size(); i++) {
      if (theta_is_pen[i]) {
        if (std::fabs(theta_value[i]) < (lambda * gamma)) {
          if (std::fabs(theta_value[i]) < DBL_EPSILON) {
            regularizer_value_i = 0;
          } else {
            regularizer_value_i = 
              lambda * (std::fabs(theta_value[i]) - std::pow(theta_value[i], 2) / (2 * lambda * gamma));
          }
        } else {
          regularizer_value_i = (std::pow(lambda, 2) * gamma) / 2;
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
    for (i = 0; i < theta_name.size(); i++) {
      if (theta_is_pen[i]) {
        if ((theta_value[i] <= (lambda * gamma)) & (theta_value[i] > 0)) {
          regularizer_gradient(i, 0) = lambda - (theta_value[i] / gamma);
        } else if ((- theta_value[i] <= (lambda * gamma)) & (theta_value[i] < 0)) {
          regularizer_gradient(i, 0) = - lambda - (theta_value[i] / gamma);
        } else if ((theta_value[i] > (lambda * gamma)) | ((- theta_value[i]) > (lambda * gamma))) {
          regularizer_gradient(i, 0) = 0;
        } else {
          regularizer_gradient(i, 0) = lambda;
        }
      } else {
        regularizer_gradient(i, 0) = 0;
      }
    }
  } else {
    regularizer_gradient = Eigen::MatrixXd::Zero(theta_name.size(), 1);
  }
}

void lslxOptimizer::update_objective_value() {
  objective_value = loss_value + regularizer_value;
}

void lslxOptimizer::update_objective_gradient() {
  int i;
  for (i = 0; i < theta_name.size(); i++) {
    if (std::fabs(theta_value[i]) > DBL_EPSILON) {
      objective_gradient(i, 0) = loss_gradient(i, 0) + regularizer_gradient(i, 0);
    } else {
      objective_gradient(i, 0) = sign(loss_gradient(i, 0)) * 
        std::fmax((std::fabs(loss_gradient(i, 0)) - lambda), 0);
    }
  }
}

void lslxOptimizer::update_theta_direction() {
  theta_direction = Rcpp::rep(0.0, theta_name.size());
  Rcpp::NumericVector z = Rcpp::rep(0.0, theta_name.size());
  Eigen::MatrixXd g = loss_gradient;
  Eigen::MatrixXd h = loss_expected_hessian;
  double z_r, z_l;
  int i, j;
  for (i = 0; i < theta_name.size(); i++) {
    h(i, i) = h(i, i) + ridge_hessian;
  }
  
  for (i = 0; i < iter_in_max; i++) {
    for (j = 0; j < theta_name.size(); j++) {
      Eigen::Map<Eigen::VectorXd> d(Rcpp::as< Eigen::Map <Eigen::VectorXd> >(theta_direction));
      double g_ij = g(j, 0) + (h * d)(j);
      double h_ij = h(j, j);
      if (lambda > DBL_EPSILON) {
        if (theta_is_free[j]) {
          z[j] = (-g_ij / h_ij);
        } else if (theta_is_pen[j]) {
          z_r = ((theta_value[j] + theta_direction[j]) / gamma - g_ij - lambda) / (h_ij - (1.0 / gamma));
          z_l = ((theta_value[j] + theta_direction[j]) / gamma - g_ij + lambda) / (h_ij - (1.0 / gamma));
          if (z_r >= - (theta_value[j] + theta_direction[j])) {
            if (z_r >= (lambda * gamma - (theta_value[j] + theta_direction[j]))) {
              z[j] = (-g_ij / h_ij);
            } else {
              z[j] = z_r;
            }
          } else if (z_l <= -(theta_value[j] + theta_direction[j])) {
            if (z_l <= (-lambda * gamma - (theta_value[j] + theta_direction[j]))) {
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
    if ((loss_expected_hessian.diagonal().array() * 
        Rcpp::as<Eigen::VectorXd>(z).array().abs()).maxCoeff() < tol_in) {
      break;
    }
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
  double regularizer_value_0;
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


void lslxOptimizer::optimize_theta_value() {
  Rcpp::NumericVector objective_gradient_abs(theta_name.size());

  update_coefficient_matrice();
  update_implied_moment();
    
  update_weight_normal();
  update_moment_gradient();
    
  update_loss_value();
  update_loss_gradient();
  update_loss_expected_hessian();
    
  update_regularizer_value();
  update_objective_value();
  
  int i, j;  
    for (i = 1; i <= iter_out_max; i++) {
      update_theta_direction();
      update_theta_value();
      update_theta_start();
      
      update_weight_normal();
      update_moment_gradient();
      
      update_loss_gradient();
      update_loss_expected_hessian();
      
      update_regularizer_gradient();
      update_objective_gradient();


      for (j = 0; j < theta_name.size(); j++) {
        if (theta_is_free[j] | theta_is_pen[j]) {
          objective_gradient_abs[j] = std::fabs(objective_gradient(j, 0));
        } else {
          objective_gradient_abs[j] = - INFINITY;
        }
      }
      objective_gradient_abs_max = Rcpp::max(objective_gradient_abs);
      n_iter_out = i;
      
      if ((objective_gradient_abs_max < tol_out) | (i == iter_out_max)) {
        update_numerical_condition();
        update_goodness_of_fit();
        break;
      }
    }
}


void lslxOptimizer::update_numerical_condition() {
  Rcpp::NumericVector objective_hessian_diagonal(theta_name.size());
  Rcpp::IntegerVector theta_is_effective(theta_name.size());
  int j;
  for (j = 0; j < theta_name.size(); j++) {
    if (theta_is_free[j]) {
      objective_hessian_diagonal[j] = loss_expected_hessian(j, j) + ridge_hessian;
      theta_is_effective[j] = 1;
    } else if (theta_is_pen[j]) {
      objective_hessian_diagonal[j] = loss_expected_hessian(j, j) + ridge_hessian - (1 / gamma);
      if (std::abs(theta_value[j]) > DBL_EPSILON) {
        theta_is_effective[j] = 1;
      } else {
        theta_is_effective[j] = 0;
      }
    } else {
      objective_hessian_diagonal[j] = INFINITY;
      theta_is_effective[j] = 0;
    }
  }
  objective_hessian_convexity = Rcpp::min(objective_hessian_diagonal);
  n_nonzero_coefficient = Rcpp::sum(theta_is_effective); 
  degree_of_freedom = n_group * n_moment - n_nonzero_coefficient;
}

void lslxOptimizer::update_goodness_of_fit() {
  aic = loss_value + (2.0 / double(total_sample_size)) * double(n_nonzero_coefficient);
  aic3 = loss_value + (3.0 / double(total_sample_size)) * double(n_nonzero_coefficient);
  caic = loss_value + ((1 + std::log(double(total_sample_size))) / double(total_sample_size)) * double(n_nonzero_coefficient);
  
  bic = loss_value + (std::log(double(total_sample_size)) / double(total_sample_size)) * double(n_nonzero_coefficient);
  abic = loss_value + (std::log((double(total_sample_size) + 2.0) / 24.0) / double(total_sample_size)) * double(n_nonzero_coefficient);
  hbic = loss_value + (std::log(double(total_sample_size) / (2.0 * 3.1415926)) / double(total_sample_size)) * double(n_nonzero_coefficient);
  
  if ((degree_of_freedom == 0) & (std::fabs(loss_value) > DBL_EPSILON)) {
    rmsea = NAN;
  } else {
    if (std::fabs(loss_value) < DBL_EPSILON) {
      rmsea = 0;
    } else {
      rmsea = std::sqrt(n_group * std::max(((loss_value / double(degree_of_freedom)) - 
        (1 / double(total_sample_size))), 0.0)); 
    }
  }
  
  double cfi_num = std::max((double(total_sample_size) * loss_value - double(degree_of_freedom)), 0.0);
  double cfi_den = std::max(std::max(double(total_sample_size) * loss_value - double(degree_of_freedom),
                                     double(total_sample_size) * loss_value_baseline - double(degree_of_freedom_baseline)), 0.0);
  if ((cfi_num < DBL_EPSILON) & (cfi_den < DBL_EPSILON)) {
    cfi = NAN;
  } else {
    if (cfi_num < DBL_EPSILON) {
      cfi = 1;
    } else {
      cfi = 1 - cfi_num / cfi_den;
    }
  }
  
  double nnfi_num = (double(total_sample_size) * loss_value_baseline) / double(degree_of_freedom_baseline) -
    (double(total_sample_size) * loss_value) / double(degree_of_freedom);
  double nnfi_den = (double(total_sample_size) * loss_value_baseline) / double(degree_of_freedom_baseline) - 1;
  nnfi = nnfi_num / nnfi_den;
  
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


Rcpp::NumericVector lslxOptimizer::extract_theta_value() {
  return Rcpp::clone(theta_value);
}

Rcpp::NumericVector lslxOptimizer::extract_numerical_condition() {
  Rcpp::NumericVector numerical_condition = 
    Rcpp::NumericVector::create(
      _["lambda"] = lambda,
      _["gamma"] = gamma,
      _["objective_value"] = objective_value,
      _["objective_gradient_abs_max"] = objective_gradient_abs_max,
      _["objective_hessian_convexity"] = objective_hessian_convexity,
      _["n_iter_out"] = n_iter_out,
      _["n_nonzero_coefficient"] = n_nonzero_coefficient,
      _["degree_of_freedom"] = degree_of_freedom);
  return Rcpp::clone(numerical_condition);
}


Rcpp::NumericVector lslxOptimizer::extract_goodness_of_fit() {
  Rcpp::NumericVector goodness_of_fit = 
    Rcpp::NumericVector::create(
      _["loss"] = loss_value,
      _["rmsea"] = rmsea,
      _["aic"] = aic,
      _["aic3"] = aic3,
      _["caic"] = caic,
      _["bic"] = bic,
      _["abic"] = abic,
      _["hbic"] = hbic,
      _["cfi"] = cfi,
      _["nnfi"] = nnfi,
      _["srmr"] = srmr);
  return Rcpp::clone(goodness_of_fit);
}


Eigen::MatrixXd lslxOptimizer::slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx) {
  Eigen::MatrixXd y(x.rows(), col_idx.size());
  int i;
  for (i = 0; i < col_idx.size(); i++) {
    y.col(i) = x.col(col_idx[i]);
  }
  return(y);
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
  return(y);
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

// [[Rcpp::export]]
void compute_regularized_path_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List numerical_condition,
    Rcpp::List goodness_of_fit,
    Rcpp::List coefficient) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  char regularizer_type = Rcpp::as<char>(control["penalty_method"]);
  Rcpp::NumericVector lambda_grid = Rcpp::as<Rcpp::NumericVector>(control["lambda_grid"]);
  Rcpp::NumericVector gamma_grid = Rcpp::as<Rcpp::NumericVector>(control["gamma_grid"]);
  
  int i, j, idx;
  for (i = 0; i < lambda_grid.size(); i++) {
    for (j = 0; j < gamma_grid.size(); j++) {
      optimizer.set_regularizer(regularizer_type, lambda_grid[i], gamma_grid[j]);
      optimizer.optimize_theta_value();
      idx = i * gamma_grid.size() + j;
      numerical_condition[idx] = optimizer.extract_numerical_condition();
      goodness_of_fit[idx] = optimizer.extract_goodness_of_fit();
      coefficient[idx] = optimizer.extract_theta_value();
    }
  }
}

// [[Rcpp::export]]
Rcpp::List compute_coefficient_matrice_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value) {
  Rcpp::List coefficient_matrice;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
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
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value) {
  Rcpp::List implied_cov;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  implied_cov = optimizer.sigma;
  return Rcpp::wrap(implied_cov);
}

// [[Rcpp::export]]
Rcpp::List compute_implied_mean_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value) {
  Rcpp::List implied_mean;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  implied_mean = optimizer.mu;
  return Rcpp::wrap(implied_mean);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_moment_gradient_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value) {
  Eigen::MatrixXd moment_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  optimizer.set_theta_value(theta_value);

  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  optimizer.update_moment_gradient();

  moment_gradient = optimizer.moment_gradient;
  return Rcpp::wrap(moment_gradient);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_expected_fisher_information_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value) {
  Eigen::MatrixXd expected_fisher_information;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  optimizer.update_weight_normal();
  optimizer.update_moment_gradient();
  optimizer.update_loss_expected_hessian();
  expected_fisher_information = 0.5 * optimizer.loss_expected_hessian;
  return Rcpp::wrap(expected_fisher_information);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_loss_gradient_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value) {
  Eigen::MatrixXd loss_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  
  optimizer.update_weight_normal();
  optimizer.update_moment_gradient();
  optimizer.update_loss_gradient();
  loss_gradient = optimizer.loss_gradient;
  return Rcpp::wrap(loss_gradient);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_regularizer_gradient_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value,
    double lambda,
    double gamma) {
  Eigen::MatrixXd regularizer_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  optimizer.set_theta_value(theta_value);
  optimizer.set_regularizer(Rcpp::as<char>(control["penalty_method"]), lambda, gamma);
  optimizer.update_regularizer_gradient();
  regularizer_gradient = optimizer.regularizer_gradient;
  return Rcpp::wrap(regularizer_gradient);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_objective_gradient_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::NumericVector theta_value,
    double lambda,
    double gamma) {
  Eigen::MatrixXd objective_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control);
  optimizer.set_theta_value(theta_value);
  optimizer.set_regularizer(Rcpp::as<char>(control["penalty_method"]), lambda, gamma);
  
  optimizer.update_coefficient_matrice();
  optimizer.update_implied_moment();
  optimizer.update_weight_normal();
  optimizer.update_moment_gradient();

  optimizer.update_loss_gradient();
  optimizer.update_regularizer_gradient();
  optimizer.update_objective_gradient();
  objective_gradient = optimizer.objective_gradient;
  return Rcpp::wrap(objective_gradient);
}



