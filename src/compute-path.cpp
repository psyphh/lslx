// Cpp class lslxOptimizer for minimizing PL criterion
// written by Po-Hsien Huang psyphh@gmail.com

#include "lslxOptimizer.h"

// [[Rcpp::depends(RcppEigen)]]
// compute solution path
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
  Rcpp::NumericVector theta_start_zero = Rcpp::clone(optimizer.theta_start);
  Rcpp::NumericVector lambda_grid = Rcpp::as<Rcpp::NumericVector>(control["lambda_grid"]);
  Rcpp::NumericVector delta_grid = Rcpp::as<Rcpp::NumericVector>(control["delta_grid"]);
  Rcpp::List numerical_condition = Rcpp::as<Rcpp::List>(fitted_result["numerical_condition"]);
  Rcpp::List information_criterion = Rcpp::as<Rcpp::List>(fitted_result["information_criterion"]);
  Rcpp::List fit_index = Rcpp::as<Rcpp::List>(fitted_result["fit_index"]);
  Rcpp::List coefficient = Rcpp::as<Rcpp::List>(fitted_result["coefficient"]);
  
  int i, j, idx;
  for (i = 0; i < lambda_grid.size(); i++) {
    if (!optimizer.warm_start) {
      optimizer.set_theta_value(theta_start_zero);
    }
    for (j = 0; j < delta_grid.size(); j++) {
      optimizer.set_regularizer(
        Rcpp::as< Rcpp::CharacterVector >(control["regularizer_type"]), 
        lambda_grid[i], delta_grid[j]);
      optimizer.complete_estimation();
      idx = i * delta_grid.size() + j;
      coefficient[idx] = optimizer.extract_coefficient();
      numerical_condition[idx] = optimizer.extract_numerical_condition();
      information_criterion[idx] = optimizer.extract_information_criterion();
      fit_index[idx] = optimizer.extract_fit_index();
    }
  }
}



// compute stepwise solution path made by forward or backward selection
// [[Rcpp::export]]
void compute_stepwise_path_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result,
    Rcpp::List fitted_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Rcpp::NumericVector theta_start_zero = Rcpp::clone(optimizer.theta_start);
  optimizer.set_regularizer(
    Rcpp::as< Rcpp::CharacterVector >(control["regularizer_type"]), 0.0, INFINITY);
  Rcpp::IntegerVector step_grid = Rcpp::as<Rcpp::IntegerVector>(control["step_grid"]);
  Rcpp::List numerical_condition = Rcpp::as<Rcpp::List>(fitted_result["numerical_condition"]);
  Rcpp::List information_criterion = Rcpp::as<Rcpp::List>(fitted_result["information_criterion"]);
  Rcpp::List fit_index = Rcpp::as<Rcpp::List>(fitted_result["fit_index"]);
  Rcpp::List coefficient = Rcpp::as<Rcpp::List>(fitted_result["coefficient"]);
  
  int i;
  for (i = 0; i < step_grid.size(); i++) {
    if (!optimizer.warm_start) {
      optimizer.set_theta_value(theta_start_zero);
    }
    if (i == 0) {
      optimizer.complete_estimation();
      coefficient[i] = optimizer.extract_coefficient();
      numerical_condition[i] = optimizer.extract_numerical_condition();
      information_criterion[i] = optimizer.extract_information_criterion();
      fit_index[i] = optimizer.extract_fit_index();
    } else {
      optimizer.complete_searching();
      coefficient[i] = optimizer.extract_coefficient();
      numerical_condition[i] = optimizer.extract_numerical_condition();
      information_criterion[i] = optimizer.extract_information_criterion();
      fit_index[i] = optimizer.extract_fit_index();
    }
  }
}


// compute solution path made by unpenalized estimation
// [[Rcpp::export]]
void compute_none_path_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result,
    Rcpp::List fitted_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Rcpp::NumericVector theta_start_zero = Rcpp::clone(optimizer.theta_start);
  optimizer.set_regularizer(
    Rcpp::as< Rcpp::CharacterVector >(control["regularizer_type"]), 0.0, INFINITY);
  Rcpp::List numerical_condition = Rcpp::as<Rcpp::List>(fitted_result["numerical_condition"]);
  Rcpp::List information_criterion = Rcpp::as<Rcpp::List>(fitted_result["information_criterion"]);
  Rcpp::List fit_index = Rcpp::as<Rcpp::List>(fitted_result["fit_index"]);
  Rcpp::List coefficient = Rcpp::as<Rcpp::List>(fitted_result["coefficient"]);
  
  optimizer.complete_estimation();
  coefficient[0] = optimizer.extract_coefficient();
  numerical_condition[0] = optimizer.extract_numerical_condition();
  information_criterion[0] = optimizer.extract_information_criterion();
  fit_index[0] = optimizer.extract_fit_index();
}

