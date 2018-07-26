// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// compute_regularized_path_cpp
void compute_regularized_path_cpp(Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result, Rcpp::List fitted_result);
RcppExport SEXP _lslx_compute_regularized_path_cpp(SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP, SEXP fitted_resultSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fitted_result(fitted_resultSEXP);
    compute_regularized_path_cpp(reduced_data, reduced_model, control, supplied_result, fitted_result);
    return R_NilValue;
END_RCPP
}
// compute_coefficient_matrix_cpp
Rcpp::List compute_coefficient_matrix_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_coefficient_matrix_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_coefficient_matrix_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_implied_cov_cpp
Rcpp::List compute_implied_cov_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_implied_cov_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_implied_cov_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_implied_mean_cpp
Rcpp::List compute_implied_mean_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_implied_mean_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_implied_mean_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_moment_jacobian_cpp
Rcpp::NumericMatrix compute_moment_jacobian_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_moment_jacobian_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_moment_jacobian_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_bfgs_hessian_cpp
Rcpp::NumericMatrix compute_bfgs_hessian_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_bfgs_hessian_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_bfgs_hessian_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_expected_fisher_cpp
Rcpp::NumericMatrix compute_expected_fisher_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_expected_fisher_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_expected_fisher_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_observed_fisher_cpp
Rcpp::NumericMatrix compute_observed_fisher_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_observed_fisher_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_observed_fisher_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_score_acov_cpp
Rcpp::NumericMatrix compute_score_acov_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_score_acov_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_score_acov_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_loss_value_cpp
double compute_loss_value_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_loss_value_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loss_value_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_loss_gradient_cpp
Rcpp::NumericMatrix compute_loss_gradient_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_loss_gradient_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loss_gradient_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_loss_gradient_direct_cpp
Rcpp::NumericMatrix compute_loss_gradient_direct_cpp(Rcpp::NumericVector theta_value, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_loss_gradient_direct_cpp(SEXP theta_valueSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loss_gradient_direct_cpp(theta_value, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_regularizer_gradient_cpp
Rcpp::NumericMatrix compute_regularizer_gradient_cpp(Rcpp::NumericVector theta_value, double lambda, double delta, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_regularizer_gradient_cpp(SEXP theta_valueSEXP, SEXP lambdaSEXP, SEXP deltaSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_regularizer_gradient_cpp(theta_value, lambda, delta, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_objective_gradient_cpp
Rcpp::NumericMatrix compute_objective_gradient_cpp(Rcpp::NumericVector theta_value, double lambda, double delta, Rcpp::List reduced_data, Rcpp::List reduced_model, Rcpp::List control, Rcpp::List supplied_result);
RcppExport SEXP _lslx_compute_objective_gradient_cpp(SEXP theta_valueSEXP, SEXP lambdaSEXP, SEXP deltaSEXP, SEXP reduced_dataSEXP, SEXP reduced_modelSEXP, SEXP controlSEXP, SEXP supplied_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta_value(theta_valueSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_data(reduced_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type reduced_model(reduced_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type supplied_result(supplied_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_objective_gradient_cpp(theta_value, lambda, delta, reduced_data, reduced_model, control, supplied_result));
    return rcpp_result_gen;
END_RCPP
}
// compute_saturated_moment_cpp
void compute_saturated_moment_cpp(Rcpp::List y_obs, Rcpp::List w, Rcpp::List m_idx, Rcpp::List saturated_mean, Rcpp::List saturated_cov, int iter_other_max, double tol_other);
RcppExport SEXP _lslx_compute_saturated_moment_cpp(SEXP y_obsSEXP, SEXP wSEXP, SEXP m_idxSEXP, SEXP saturated_meanSEXP, SEXP saturated_covSEXP, SEXP iter_other_maxSEXP, SEXP tol_otherSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type y_obs(y_obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type w(wSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type m_idx(m_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type saturated_mean(saturated_meanSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type saturated_cov(saturated_covSEXP);
    Rcpp::traits::input_parameter< int >::type iter_other_max(iter_other_maxSEXP);
    Rcpp::traits::input_parameter< double >::type tol_other(tol_otherSEXP);
    compute_saturated_moment_cpp(y_obs, w, m_idx, saturated_mean, saturated_cov, iter_other_max, tol_other);
    return R_NilValue;
END_RCPP
}
// compute_saturated_moment_acov_response_cpp
void compute_saturated_moment_acov_response_cpp(Rcpp::List y_obs, Rcpp::List w, Rcpp::List m_idx, Rcpp::List m2_idx, Rcpp::List saturated_mean, Rcpp::List saturated_cov, Rcpp::List saturated_moment_acov);
RcppExport SEXP _lslx_compute_saturated_moment_acov_response_cpp(SEXP y_obsSEXP, SEXP wSEXP, SEXP m_idxSEXP, SEXP m2_idxSEXP, SEXP saturated_meanSEXP, SEXP saturated_covSEXP, SEXP saturated_moment_acovSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type y_obs(y_obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type w(wSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type m_idx(m_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type m2_idx(m2_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type saturated_mean(saturated_meanSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type saturated_cov(saturated_covSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type saturated_moment_acov(saturated_moment_acovSEXP);
    compute_saturated_moment_acov_response_cpp(y_obs, w, m_idx, m2_idx, saturated_mean, saturated_cov, saturated_moment_acov);
    return R_NilValue;
END_RCPP
}
// compute_saturated_moment_acov_moment_cpp
void compute_saturated_moment_acov_moment_cpp(int n_observation, Rcpp::List sample_proportion, Rcpp::List saturated_cov, Rcpp::List saturated_moment_acov);
RcppExport SEXP _lslx_compute_saturated_moment_acov_moment_cpp(SEXP n_observationSEXP, SEXP sample_proportionSEXP, SEXP saturated_covSEXP, SEXP saturated_moment_acovSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_observation(n_observationSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sample_proportion(sample_proportionSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type saturated_cov(saturated_covSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type saturated_moment_acov(saturated_moment_acovSEXP);
    compute_saturated_moment_acov_moment_cpp(n_observation, sample_proportion, saturated_cov, saturated_moment_acov);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lslx_compute_regularized_path_cpp", (DL_FUNC) &_lslx_compute_regularized_path_cpp, 5},
    {"_lslx_compute_coefficient_matrix_cpp", (DL_FUNC) &_lslx_compute_coefficient_matrix_cpp, 5},
    {"_lslx_compute_implied_cov_cpp", (DL_FUNC) &_lslx_compute_implied_cov_cpp, 5},
    {"_lslx_compute_implied_mean_cpp", (DL_FUNC) &_lslx_compute_implied_mean_cpp, 5},
    {"_lslx_compute_moment_jacobian_cpp", (DL_FUNC) &_lslx_compute_moment_jacobian_cpp, 5},
    {"_lslx_compute_bfgs_hessian_cpp", (DL_FUNC) &_lslx_compute_bfgs_hessian_cpp, 5},
    {"_lslx_compute_expected_fisher_cpp", (DL_FUNC) &_lslx_compute_expected_fisher_cpp, 5},
    {"_lslx_compute_observed_fisher_cpp", (DL_FUNC) &_lslx_compute_observed_fisher_cpp, 5},
    {"_lslx_compute_score_acov_cpp", (DL_FUNC) &_lslx_compute_score_acov_cpp, 5},
    {"_lslx_compute_loss_value_cpp", (DL_FUNC) &_lslx_compute_loss_value_cpp, 5},
    {"_lslx_compute_loss_gradient_cpp", (DL_FUNC) &_lslx_compute_loss_gradient_cpp, 5},
    {"_lslx_compute_loss_gradient_direct_cpp", (DL_FUNC) &_lslx_compute_loss_gradient_direct_cpp, 5},
    {"_lslx_compute_regularizer_gradient_cpp", (DL_FUNC) &_lslx_compute_regularizer_gradient_cpp, 7},
    {"_lslx_compute_objective_gradient_cpp", (DL_FUNC) &_lslx_compute_objective_gradient_cpp, 7},
    {"_lslx_compute_saturated_moment_cpp", (DL_FUNC) &_lslx_compute_saturated_moment_cpp, 7},
    {"_lslx_compute_saturated_moment_acov_response_cpp", (DL_FUNC) &_lslx_compute_saturated_moment_acov_response_cpp, 7},
    {"_lslx_compute_saturated_moment_acov_moment_cpp", (DL_FUNC) &_lslx_compute_saturated_moment_acov_moment_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_lslx(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
