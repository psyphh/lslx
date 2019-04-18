## \code{$extract_specification()} returns a \code{data.frame} of model specification. ##
lslx$set("public",
         "extract_specification",
         function() {
           specification <-
             private$model$specification
           return(specification)
         })

## \code{$extract_saturated_cov()} returns a \code{list} of saturated sample covariance matrice(s). ##
lslx$set("public",
         "extract_saturated_cov",
         function() {
           saturated_cov <-
             private$fitting$reduced_data$saturated_cov
           return(saturated_cov)
         })

## \code{$extract_saturated_mean()} returns a \code{list} of saturated sample mean vector(s). ##
lslx$set("public",
         "extract_saturated_mean",
         function() {
           saturated_mean <-
             private$fitting$reduced_data$saturated_mean
           return(saturated_mean)
         })

## \code{$extract_saturated_moment_acov()} returns a \code{list} of asymptotic covariance matrice(s) of saturated moments. ##
lslx$set("public",
         "extract_saturated_moment_acov",
         function() {
           saturated_moment_acov <-
             private$fitting$reduced_data$saturated_moment_acov
           return(saturated_moment_acov)
         })

## \code{$extract_lambda_grid()} returns a \code{numeric} of lambda grid. ##
lslx$set("public",
         "extract_lambda_grid",
         function() {
           lambda_grid <-
             private$fitting$control$lambda_grid
           return(lambda_grid)
         })

## \code{$extract_delta_grid()} returns a \code{numeric} of delta grid. ##
lslx$set("public",
         "extract_delta_grid",
         function() {
           delta_grid <-
             private$fitting$control$delta_grid
           return(delta_grid)
         })


## \code{$extract_weight_matrix()} returns a \code{list} of weight matrix. ##
lslx$set("public",
         "extract_weight_matrix",
         function() {
           if (!(private$fitting$control$loss %in% c("uls", "dwls", "wls"))) {
             stop("Weight matrix can be only extracted when 'loss' = 'uls', 'dwls', or 'wls'.")
           } else {
             weight_matrix <-
               private$fitting$control$weight_matrix
           }
           return(weight_matrix)
         })

## \code{$extract_penalty_level()} returns a \code{character} of the index name of the optimal penalty level. ##
lslx$set("public",
         "extract_penalty_level",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           if (!include_faulty) {
             idx_convergent <-
               which(private$fitting$fitted_result$is_convergent)
             if (length(idx_convergent) == 0) {
               stop(
                 "PL estimate under EACH penalty level is derived under nonconverged result. \n",
                 "  To include faulty results, please set 'include_faulty = TRUE'. \n"
               )
             }
             idx_convex <-
               which(private$fitting$fitted_result$is_convex)
             if (length(idx_convex) == 0) {
               stop(
                 "PL estimate under EACH penalty level is derived under nonconvex approximated hessian. \n",
                 "  To include faulty results, please set 'include_faulty = TRUE'. \n"
               )
             }
           } else {
             idx_convergent <-
               seq_len(length(private$fitting$fitted_result$numerical_condition))
             idx_convex <-
               seq_len(length(private$fitting$fitted_result$numerical_condition))
           }
           idx_used <-
             intersect(x = idx_convergent, y = idx_convex)
           if (length(idx_used) == 0) {
             stop(
               "PL estimate under EACH penalty/convex level is derive under problematic optimization.\n",
               "  To include faulty results, please set 'include_faulty = TRUE'. \n"
             )
           }
           
           if (length(private$fitting$fitted_result$numerical_condition) == 1) {
             penalty_level <-
               names(private$fitting$fitted_result$numerical_condition[idx_used])
           } else {
             
             if (private$fitting$control$regularizer) {
               if (missing(selector) & missing(lambda) & missing(delta)) {
                 if (length(private$fitting$fitted_result$numerical_condition) > 1) {
                   stop(
                     "Argument 'selector', 'lambda', and 'delta' cannot be all empty if there are many regularization levels."
                   )
                 }
               }
             } else {}
             
             if (private$fitting$control$searcher) {
               if (missing(selector) & missing(step)) {
                 if (length(private$fitting$fitted_result$numerical_condition) > 1) {
                   stop(
                     "Argument 'selector' and 'step' cannot be all empty if there are many searching steps."
                   )
                 }
               }
             } else {}
             
             if (!missing(selector)) {
               if (length(selector) > 1) {
                 stop("The length of argument 'selector' can be only one.\n")
               }
               if (!(selector %in% names(private$fitting$fitted_result$information_criterion[[1]]))) {
                 stop(
                   "Argument 'selector' is unrecognized.\n",
                   "  Selector currently recognized by 'lslx' is \n  ",
                   do.call(paste, as.list(
                     names(
                       private$fitting$fitted_result$information_criterion[[1]]
                     )
                   )),
                   "."
                 )
               }
               if ((selector %in% c("raic", "raic3", "rcaic", "rbic", "rabic", "rhbic")) &
                   (!private$fitting$control$response)) {
                 stop(
                   "When lslx object is initialized via moments,",
                   " 'raic', 'raic3', 'rcaic', 'rbic', 'rabic', and 'rhbic' are not available."
                 )
               }
               
               if (private$fitting$control$regularizer) {
                 if ((!missing(lambda)) | (!missing(delta))) {
                   stop("When 'selector' is specified, 'lambda' or 'delta' will not be used.\n")
                 }
               } else {}
               if (private$fitting$control$searcher) {
                 if ((!missing(step))) {
                   stop("When 'selector' is specified, 'step' will not be used.\n")
                 }
               } else {}
               
               penalty_level <-
                 sapply(
                   X = selector,
                   FUN = function(selector_i) {
                     information_criterion_i <- sapply(
                       X = private$fitting$fitted_result$information_criterion,
                       FUN = function(information_criterion_j) {
                         getElement(object = information_criterion_j,
                                    name = selector_i)
                       }
                     )
                     penalty_level_i <-
                       names(which.min(information_criterion_i[idx_used]))
                     return(penalty_level_i)
                   }
                 )
             } else {
               if (private$fitting$control$regularizer) {
                 if (missing(lambda)) {
                   stop("When 'selector' is not specified, 'lambda' cannot be empty.")
                 }
                 if (missing(delta)) {
                   if (private$fitting$control$penalty_method %in% c("elastic_net", "mcp")) {
                     stop(
                       "When 'selector' is not specified, 'delta' cannot be empty."
                     )
                   } else if (private$fitting$control$penalty_method == "lasso") {
                     delta <- 1
                   } else if (private$fitting$control$penalty_method == "ridge") {
                     delta <- 0
                   } else{}
                 }
                 penalty_used_split <-
                   strsplit(
                     x = names(private$fitting$fitted_result$numerical_condition[idx_used]),
                     split = "=|/"
                   )
                 penalty_used <-
                   sapply(
                     X = penalty_used_split,
                     FUN = function(penalty_used_split_i) {
                       penalty_used_i <-
                         as.numeric(c(penalty_used_split_i[2],
                                      penalty_used_split_i[4]))
                       return(penalty_used_i)
                     }
                   )
                 delta_used <- penalty_used[2, ]
                 if (delta %in% delta_used) {
                   delta_nearest <- delta
                 } else {
                   if (delta == Inf) {
                     delta_nearest <- max(delta_used)
                   } else {
                     delta_nearest <-
                       max(abs(unique(delta_used)[which.min(abs(unique(delta_used) - delta))]))
                   }
                 }
                 lambda_used <-
                   penalty_used[1, penalty_used[2, ] == delta_nearest]
                 if (lambda %in% lambda_used) {
                   lambda_nearest <- lambda
                 } else {
                   if (lambda == Inf) {
                     lambda_nearest <- max(lambda_used)
                   } else {
                     lambda_nearest <-
                       max(abs(unique(lambda_used)[which.min(abs(unique(lambda_used) - lambda))]))
                   }
                 }
                 penalty_level <-
                   paste0("ld=", lambda_nearest, "/", "gm=", delta_nearest)
               } else {}
               if (private$fitting$control$searcher) {
                 if (missing(step)) {
                   stop("When 'selector' is not specified, 'step' cannot be empty.")
                 }
                 if (private$fitting$control$penalty_method == "forward") {
                   step_nearest <-
                     min(abs(unique(private$fitting$control$step_grid)[
                       which.min(abs(unique(private$fitting$control$step_grid) - step))]))
                 } else if (private$fitting$control$penalty_method == "backward") {
                   step_nearest <-
                     max(abs(unique(private$fitting$control$step_grid)[
                       which.min(abs(unique(private$fitting$control$step_grid) - step))]))
                 } else {}
                 penalty_level <- paste0("step=", step_nearest)
               } else {}
             }
           }
           return(penalty_level)
         })

## \code{$extract_coefficient_indicator()} returns a \code{logic} by testing the types of coefficients. ##
lslx$set("public",
         "extract_coefficient_indicator",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           if (!(type %in% c(
             "default",
             "all",
             "fixed",
             "free",
             "pen",
             "effective",
             "selected"
           ))) {
             stop(
               "Argument 'type' can be only either 'default', 'all', 'fixed', 'free', 'pen', 'active', or 'selected'."
             )
           }
           if (type == "default") {
             type <- "all"
           }
           if (type == "all") {
             coefficient_indicator <- rep(T, length(coefficient))
           } else if (type == "fixed") {
             coefficient_indicator <-
               !(
                 private$fitting$reduced_model$theta_is_free |
                   private$fitting$reduced_model$theta_is_pen
               )
           } else if (type == "free") {
             coefficient_indicator <-
               private$fitting$reduced_model$theta_is_free
           } else if (type == "pen") {
             coefficient_indicator <-
               private$fitting$reduced_model$theta_is_pen
           } else if (type == "effective") {
             coefficient_indicator <-
               private$fitting$reduced_model$theta_is_free |
               (private$fitting$reduced_model$theta_is_pen &
                  coefficient != 0)
           } else if (type == "selected") {
             coefficient_indicator <-
               private$fitting$reduced_model$theta_is_pen &
               coefficient != 0
           } else {}
           return(coefficient_indicator)
         })


## \code{$extract_numerical_condition()} returns a \code{numeric} of the numerical conditions. ##
lslx$set("public",
         "extract_numerical_condition",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           numerical_condition <-
             private$fitting$fitted_result$numerical_condition[[penalty_level]]
           return(numerical_condition)
         })

## \code{$extract_information_criterion()} returns a \code{numeric} of the values of information criteria. ##
lslx$set("public",
         "extract_information_criterion",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           information_criterion <-
             private$fitting$fitted_result$information_criterion[[penalty_level]]
           return(information_criterion)
         })

## \code{$extract_fit_index()} returns a \code{numeric} of the values of fit indices. ##
lslx$set("public",
         "extract_fit_index",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           fit_index <-
             private$fitting$fitted_result$fit_index[[penalty_level]]
           return(fit_index)
         })


## \code{$extract_cv_error()} returns a \code{numeric} of the values of cv errors. ##
lslx$set("public",
         "extract_cv_error",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           cv_error <-
             private$fitting$fitted_result$cv_error[[penalty_level]]
           return(cv_error)
         })


## \code{$extract_coefficient()} returns a \code{numeric} of estimates of the coefficients. ##
lslx$set("public",
         "extract_coefficient",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             coefficient <- coefficient[coefficient_indicator]
           }
           return(coefficient)
         })



## \code{$extract_debiased_coefficient()} returns a \code{numeric} of debiased estimates of the coefficients. ##
lslx$set("public",
         "extract_debiased_coefficient",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           coefficient <-
             self$extract_coefficient(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           debiased_coefficient <- coefficient
           if (private$fitting$control$regularizer) {
             is_effective <-
               self$extract_coefficient_indicator(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 type = "effective",
                 include_faulty = include_faulty
               )
             is_selected <-
               self$extract_coefficient_indicator(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 type = "selected",
                 include_faulty = include_faulty
               )
             if (any(is_selected)) {
               penalty_level <-
                 self$extract_penalty_level(
                   selector = selector,
                   lambda = lambda,
                   delta = delta,
                   step = step,
                   include_faulty = include_faulty
                 )
               lambda_ <- as.numeric(strsplit(x = penalty_level,
                                              split = "=|/")[[1]][2])
               delta_ <- as.numeric(strsplit(x = penalty_level,
                                             split = "=|/")[[1]][4])
               regularizer_gradient <-
                 compute_regularizer_gradient_cpp(
                   theta_value = coefficient,
                   lambda = lambda_,
                   delta = delta_,
                   reduced_data = private$fitting$reduced_data,
                   reduced_model = private$fitting$reduced_model,
                   control = private$fitting$control,
                   supplied_result = private$fitting$supplied_result
                 )
               observed_information <-
                 2 * self$extract_observed_information(
                   selector = selector,
                   lambda = lambda,
                   delta = delta,
                   step = step,
                   include_faulty = include_faulty
                 )
               observed_information_inv <-
                 matrix(0, length(coefficient), length(coefficient))
               observed_information_inv[is_effective, is_effective] <-
                 solve(observed_information[is_effective, is_effective])
               debiased_coefficient[is_effective] <-
                 coefficient[is_effective] +
                 observed_information_inv[is_effective, is_effective, drop = FALSE] %*%
                 (regularizer_gradient[is_effective, 1, drop = FALSE])
             }
           } else {}
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             debiased_coefficient <- debiased_coefficient[coefficient_indicator]
           }
           return(debiased_coefficient)
         })


## \code{$extract_implied_cov()} returns a \code{list} of model-implied covariance matrice(s). ##
lslx$set("public",
         "extract_implied_cov",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           implied_cov <-
             compute_implied_cov_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           implied_cov <-
             lapply(
               X = implied_cov,
               FUN = function(implied_cov_i) {
                 rownames(implied_cov_i) <- private$model$name_response
                 colnames(implied_cov_i) <-
                   private$model$name_response
                 return(implied_cov_i)
               }
             )
           names(implied_cov) <- private$model$level_group
           return(implied_cov)
         })

## \code{$extract_implied_mean()} returns a \code{list} of model-implied mean vector(s). ##
lslx$set("public",
         "extract_implied_mean",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           implied_mean <-
             compute_implied_mean_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           implied_mean <-
             lapply(
               X = implied_mean,
               FUN = function(implied_mean_i) {
                 rownames(implied_mean_i) <- private$model$name_response
                 return(implied_mean_i)
               }
             )
           names(implied_mean) <- private$model$level_group
           return(implied_mean)
         })

## \code{$extract_residual_cov()} returns a \code{list} of residual matrice(s) of covariance. ##
lslx$set("public",
         "extract_residual_cov",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           implied_cov <-
             self$extract_implied_cov(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           residual_cov <-
             mapply(
               FUN = function(implied_cov_i,
                              saturated_cov_i) {
                 residual_cov_i <-
                   implied_cov_i - saturated_cov_i
                 return(residual_cov_i)
               },
               implied_cov,
               private$fitting$reduced_data$saturated_cov,
               SIMPLIFY = FALSE,
               USE.NAMES = TRUE
             )
           return(residual_cov)
         })

## \code{$extract_residual_mean()} returns a \code{list} of residual vector(s) of mean. ##
lslx$set("public",
         "extract_residual_mean",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           implied_mean <-
             self$extract_implied_mean(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           residual_mean <-
             mapply(
               FUN = function(implied_mean_i,
                              saturated_mean_i) {
                 residual_mean_i <-
                   implied_mean_i - saturated_mean_i
                 return(residual_mean_i)
               },
               implied_mean,
               private$fitting$reduced_data$saturated_mean,
               SIMPLIFY = FALSE,
               USE.NAMES = TRUE
             )
           return(residual_mean)
         })

## \code{$extract_coefficient_matrix()} returns a \code{list} of coefficient matrice(s) specified by \code{block}. ##
lslx$set("public",
         "extract_coefficient_matrix",
         function(selector,
                  lambda,
                  delta,
                  step,
                  block,
                  include_faulty = FALSE) {
           if (missing(block)) {
             stop("Argument 'block' is missing.")
           }
           if (length(block) > 1) {
             stop("The length of argument 'block' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           coefficient_matrix <-
             compute_coefficient_matrix_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           
           if (block %in% c("f<-1", "y<-1")) {
             selected_matrix <- coefficient_matrix$alpha
             col_names <- "1"
           } else {
             if (block %in% c("f<-f", "f<-y", "y<-f", "y<-y")) {
               selected_matrix <- coefficient_matrix$beta
             } else if (block %in% c("f<->f", "f<->y", "y<->f", "y<->y")) {
               selected_matrix <- coefficient_matrix$phi
             } else {
               stop(
                 "Argument 'block' is unrecognized. It must be one of the following:\n",
                 "  'f<-1', 'y<-1', 'f<-f', 'f<-y', 'y<-f', 'f<-f', 'f<->f', 'f<->y', 'y<->f', 'f<->f'."
               )
             }
             col_names <- private$model$name_eta
           }
           row_select <-
             strsplit(block, split = "<-|<->|->")[[1]][1]
           col_select <-
             strsplit(block, split = "<-|<->|->")[[1]][2]
           if (row_select == "f") {
             row_select <- private$model$name_factor
           } else if (row_select == "y") {
             row_select <- private$model$name_response
           }
           if (col_select == "f") {
             col_select <- private$model$name_factor
           } else if (col_select == "y") {
             col_select <- private$model$name_response
           }
           selected_matrix <-
             lapply(
               X = selected_matrix,
               FUN = function(selected_matrix_i) {
                 rownames(selected_matrix_i) <- private$model$name_eta
                 colnames(selected_matrix_i) <- col_names
                 return(selected_matrix_i)
               }
             )
           coefficient_matrix_block <-
             lapply(
               X = selected_matrix,
               FUN = function(selected_matrix_i) {
                 coefficient_matrix_block_i <-
                   selected_matrix_i[row_select,
                                     col_select,
                                     drop = FALSE]
                 return(coefficient_matrix_block_i)
               }
             )
           names(coefficient_matrix_block) <-
             private$model$level_group
           return(coefficient_matrix_block)
         })


## \code{$extract_moment_jacobian()} returns a \code{matrix} of Jacobian of moment structure. ##
lslx$set("public",
         "extract_moment_jacobian",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           moment_jacobian <-
             compute_model_jacobian_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(moment_jacobian) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             moment_jacobian <-
               moment_jacobian[, coefficient_indicator, drop = FALSE]
           }
           return(moment_jacobian)
         })

## \code{$extract_expected_information()} returns a \code{matrix} of the expected Fisher information matrix. ##
lslx$set("public",
         "extract_expected_information",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           expected_information <-
             compute_expected_information_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(expected_information) <-
             rownames(private$model$specification)
           rownames(expected_information) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             expected_information <- expected_information[coefficient_indicator,
                                                          coefficient_indicator,
                                                          drop = FALSE]
           }
           return(expected_information)
         })


## \code{$extract_observed_information()} returns a \code{matrix} of the observed Fisher information matrix. ##
lslx$set("public",
         "extract_observed_information",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           observed_information <-
             compute_observed_information_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(observed_information) <-
             rownames(private$model$specification)
           rownames(observed_information) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             observed_information <- observed_information[coefficient_indicator,
                                                          coefficient_indicator,
                                                          drop = FALSE]
           }
           return(observed_information)
         })

## \code{$extract_bfgs_hessian()} returns a \code{matrix} of the BFGS Hessian matrix. ##
lslx$set("public",
         "extract_bfgs_hessian",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           bfgs_hessian <-
             compute_bfgs_hessian_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(bfgs_hessian) <-
             rownames(private$model$specification)
           rownames(bfgs_hessian) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             bfgs_hessian <- bfgs_hessian[coefficient_indicator,
                                          coefficient_indicator,
                                          drop = FALSE]
           }
           return(bfgs_hessian)
         })

## \code{$extract_score_acov()} returns a \code{matrix} of the asymptotic covariance of scores. ##
lslx$set("public",
         "extract_score_acov",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           score_acov <-
             compute_score_acov_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(score_acov) <-
             rownames(private$model$specification)
           rownames(score_acov) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             score_acov <- score_acov[coefficient_indicator,
                                      coefficient_indicator,
                                      drop = FALSE]
           }
           return(score_acov)
         })

## \code{$extract_coefficient_acov()} returns a \code{matrix} of the asymptotic covariance of coefficients. ##
lslx$set("public",
         "extract_coefficient_acov",
         function(selector,
                  lambda,
                  delta,
                  step,
                  standard_error = "default",
                  type = "default",
                  include_faulty = FALSE) {
           if (!(
             standard_error %in% c(
               "default",
               "sandwich",
               "observed_information",
               "expected_information"
             )
           )) {
             stop(
               "Argument 'standard_error' can be only either 'default', 'sandwich', 'observed_information', or 'expected_information'."
             )
           }
           if (standard_error == "default") {
             if (private$fitting$control$response) {
               standard_error <- "sandwich"
             } else {
               standard_error <- "observed_information"
             }
           }
           coefficient <-
             self$extract_coefficient(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           is_effective <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = "effective",
               include_faulty = include_faulty
             )
           coefficient_acov <-
             matrix(NA, length(coefficient), length(coefficient))
           if (standard_error == "sandwich") {
             score_acov <-
               self$extract_score_acov(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 include_faulty = include_faulty
               )
             observed_information <-
               self$extract_observed_information(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 include_faulty = include_faulty
               )
             observed_information_pinv <-
               solve(observed_information[is_effective, is_effective])
             coefficient_acov[is_effective, is_effective] <-
               (observed_information_pinv %*%
                  score_acov[is_effective, is_effective] %*%
                  observed_information_pinv)
           } else if (standard_error == "expected_information") {
             expected_information <-
               self$extract_expected_information(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 include_faulty = include_faulty
               )
             coefficient_acov[is_effective, is_effective] <-
               solve(expected_information[is_effective, is_effective]) /
               private$fitting$reduced_data$n_observation
           } else if (standard_error == "observed_information") {
             observed_information <-
               self$extract_observed_information(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 include_faulty = include_faulty
               )
             coefficient_acov[is_effective, is_effective] <-
               solve(observed_information[is_effective, is_effective]) /
               private$fitting$reduced_data$n_observation
           } else {}
           colnames(coefficient_acov) <-
             rownames(private$model$specification)
           rownames(coefficient_acov) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             coefficient_acov <- coefficient_acov[coefficient_indicator,
                                                  coefficient_indicator,
                                                  drop = FALSE]
           }
           attr(coefficient_acov, "standard_error") <- standard_error
           return(coefficient_acov)
         })

## \code{$extract_loss_gradient()} returns a \code{matrix} of the gradient of loss function. ##
lslx$set("public",
         "extract_loss_gradient",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           loss_gradient <-
             compute_loss_gradient_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           rownames(loss_gradient) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             loss_gradient <- loss_gradient[coefficient_indicator, ,
                                            drop = FALSE]
           }
           return(loss_gradient)
         })

## \code{$extract_regularizer_gradient()} returns a \code{matrix} of the sub-gradient of regularizer. ##
lslx$set("public",
         "extract_regularizer_gradient",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           if (private$fitting$control$regularizer) {
             penalty_level <-
               self$extract_penalty_level(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 include_faulty = include_faulty
               )
             coefficient <-
               private$fitting$fitted_result$coefficient[[penalty_level]]
             lambda_ <- as.numeric(strsplit(x = penalty_level,
                                            split = "=|/")[[1]][2])
             delta_ <- as.numeric(strsplit(x = penalty_level,
                                           split = "=|/")[[1]][4])
             regularizer_gradient <-
               compute_regularizer_gradient_cpp(
                 theta_value = coefficient,
                 lambda = lambda_,
                 delta = delta_,
                 reduced_data = private$fitting$reduced_data,
                 reduced_model = private$fitting$reduced_model,
                 control = private$fitting$control,
                 supplied_result = private$fitting$supplied_result
               )
           } else {
             regularizer_gradient <- 
               matrix(0, nrow = private$fitting$reduced_model$n_theta, ncol = 1)
           }
           rownames(regularizer_gradient) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             regularizer_gradient <-
               regularizer_gradient[coefficient_indicator, ,
                                    drop = FALSE]
           }
           return(regularizer_gradient)
         })

## \code{$extract_objective_gradient()} returns a \code{matrix} of the sub-gradient of objective function. ##
lslx$set("public",
         "extract_objective_gradient",
         function(selector,
                  lambda,
                  delta,
                  step,
                  type = "default",
                  include_faulty = FALSE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               include_faulty = include_faulty
             )
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           if (private$fitting$control$regularizer) {
             lambda_ <- as.numeric(strsplit(x = penalty_level,
                                            split = "=|/")[[1]][2])
             delta_ <- as.numeric(strsplit(x = penalty_level,
                                           split = "=|/")[[1]][4])
             objective_gradient <-
               compute_objective_gradient_cpp(
                 theta_value = coefficient,
                 lambda = lambda_,
                 delta = delta_,
                 reduced_data = private$fitting$reduced_data,
                 reduced_model = private$fitting$reduced_model,
                 control = private$fitting$control,
                 supplied_result = private$fitting$supplied_result
               )
           } else {
             objective_gradient <-
               compute_loss_gradient_cpp(
                 theta_value = coefficient,
                 reduced_data = private$fitting$reduced_data,
                 reduced_model = private$fitting$reduced_model,
                 control = private$fitting$control,
                 supplied_result = private$fitting$supplied_result
               )
           } 
           rownames(objective_gradient) <-
             rownames(private$model$specification)
           coefficient_indicator <-
             self$extract_coefficient_indicator(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               type = type,
               include_faulty = include_faulty
             )
           if (!(all(coefficient_indicator))) {
             objective_gradient <- 
               objective_gradient[coefficient_indicator, ,
                                  drop = FALSE]
           }
           return(objective_gradient)
         })
