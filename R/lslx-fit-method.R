## \code{$fit()} fits the specified model to data by minimizing a penalized ML loss function. ##
lslx$set("public",
         "fit",
         function(penalty_method = "none",
                  lambda_grid = "default",
                  delta_grid = "default",
                  algorithm = "default",
                  missing_method = "default",
                  start_method = "default",
                  subset = NULL,
                  cv_fold = 1L,
                  iter_out_max = 100L,
                  iter_in_max = 30L,
                  iter_other_max = 500L,
                  iter_armijo_max = 100L,
                  tol_out = 1e-3,
                  tol_in = 1e-3,
                  tol_other = 1e-7,
                  step_size = 0.5,
                  armijo = 1e-5,
                  ridge_cov = 0,
                  ridge_hessian = 1e-4,
                  positive_diag = TRUE,
                  verbose = TRUE) {
           if (!(penalty_method %in% c("none", "lasso", "mcp"))) {
             stop("Argument 'penalty_method' can be only either 'none', 'lasso', or 'mcp'.")
           }
           if (!is.numeric(lambda_grid)) {
             if (!is.character(lambda_grid)) {
               stop("Argument 'lambda_grid' can be only a numeric vector or set as 'default'.")
             } else if (is.character(lambda_grid) &
                        (length(lambda_grid) != 1)) {
               stop("Argument 'lambda_grid' can be only a numeric vector or set as 'default'.")
             } else if (is.character(lambda_grid) &
                        (length(lambda_grid) == 1)) {
               if (lambda_grid != "default") {
                 stop("Argument 'lambda_grid' can be only a numeric vector or set as 'default'.")
               }
             }
           }
           if (!is.numeric(delta_grid)) {
             if (!is.character(delta_grid)) {
               stop("Argument 'delta_grid' can be only a numeric vector or set as 'default'.")
             } else if (is.character(delta_grid) &
                        (length(delta_grid) != 1)) {
               stop("Argument 'delta_grid' can be only a numeric vector or set as 'default'.")
             } else if (is.character(delta_grid) &
                        (length(delta_grid) == 1)) {
               if (delta_grid != "default") {
                 stop("Argument 'delta_grid' can be only a numeric vector or set as 'default'.")
               }
             }
           }
           if (!(algorithm %in% c("default", "bfgs", "fisher"))) {
             stop("Argument 'algorithm' can be only 'default', 'bfgs', or 'fisher'.")
           }
           if (!(missing_method %in% c("default", "two_stage", "listwise_deletion"))) {
             stop(
               "Argument 'start_method' can be only 'default', 'two_stage', or 'listwise_deletion'."
             )
           }
           if (!(start_method %in% c("default", "none", "mh", "heuristic"))) {
             stop("Argument 'start_method' can be only 'default', 'none', 'mh', or 'heuristic'.")
           }
           if (!is.null(subset)) {
             if (!(is.integer(subset) | is.logical(subset))) {
               stop("Argument 'subset' must be a integer or logical vector.")
             }
           }
           if (!(is.numeric(cv_fold) & (cv_fold > 0))) {
             stop("Argument 'cv_fold' must be a positive integer.")
           }
           if (!(is.numeric(iter_out_max) &
                 (length(iter_out_max) = 1))) {
             stop("Argument 'iter_out_max' must be a numeric vector with length one.")
           }
           if (!(is.numeric(iter_in_max) &
                 (length(iter_in_max) = 1))) {
             stop("Argument 'iter_in_max' must be a numeric vector with length one.")
           }
           if (!(is.numeric(iter_armijo_max) &
                 (length(iter_armijo_max) = 1))) {
             stop("Argument 'iter_armijo_max' must be a numeric vector with length one.")
           }
           if (!(is.numeric(tol_out) & (length(tol_out) = 1))) {
             stop("Argument 'tol_out' must be a numeric vector with length one.")
           }
           if (!(is.numeric(tol_in) & (length(tol_in) = 1))) {
             stop("Argument 'tol_in' must be a numeric vector with length one.")
           }
           if (!(is.numeric(step_size) & (length(step_size) = 1))) {
             stop("Argument 'step_size' must be a numeric vector with length one.")
           }
           if (!(is.numeric(armijo) & (length(armijo) = 1))) {
             stop("Argument 'armijo' must be a numeric vector with length one.")
           }
           if (!(is.numeric(ridge_cov) & (length(ridge_cov) = 1))) {
             stop("Argument 'ridge_cov' must be a numeric vector with length one.")
           }
           if (!(is.numeric(ridge_hessian) &
                 (length(ridge_hessian) = 1))) {
             stop("Argument 'ridge_hessian' must be a numeric vector with length one.")
           }
           if (!(is.logical(positive_diag) &
                 (length(positive_diag) = 1))) {
             stop("Argument 'positive_diag' must be a logical vector with length one.")
           }
           
           control <-
             list(
               penalty_method = penalty_method,
               lambda_grid = lambda_grid,
               delta_grid = delta_grid,
               algorithm = algorithm,
               missing_method = missing_method,
               start_method = start_method,
               subset = subset,
               cv_fold = cv_fold,
               iter_out_max = iter_out_max,
               iter_in_max = iter_in_max,
               iter_other_max = iter_other_max,
               iter_armijo_max = iter_armijo_max,
               tol_out = tol_out,
               tol_in = tol_in,
               tol_other = tol_other,
               step_size = step_size,
               armijo = armijo,
               ridge_cov = ridge_cov,
               ridge_hessian = ridge_hessian,
               positive_diag = positive_diag
             )
           private$fitting <-
             lslxFitting$new(model = private$model,
                             data = private$data,
                             control = control)
           compute_regularized_path_cpp(
             private$fitting$reduced_data,
             private$fitting$reduced_model,
             private$fitting$control,
             private$fitting$supplied_result,
             private$fitting$fitted_result
           )
           
           name_grid <-
             paste0(
               "ld=",
               sapply(
                 X = private$fitting$fitted_result$numerical_condition,
                 FUN = function(x) {
                   getElement(x, "lambda")
                 }
               ),
               "/",
               "gm=",
               sapply(
                 X = private$fitting$fitted_result$numerical_condition,
                 FUN = function(x) {
                   getElement(x, "delta")
                 }
               )
             )
           
           if (private$fitting$control$cv_fold > 1L) {
             if (private$fitting$control$response) {
               control <- private$fitting$control
               control$start_method <- "none"       
               control$cv_fold <- 1L
               cv_error <-
                 sapply(
                   X = 1:private$fitting$control$cv_fold,
                   FUN = function(i) {
                     subset_train_i <- 
                       private$fitting$control$subset[private$fitting$control$cv_idx != i]
                     control$subset <- subset_train_i
                     fitting_train_i <-
                       lslxFitting$new(model = private$model,
                                       data = private$data,
                                       control = control)
                     fitting_train_i$supplied_result$fitted_start <- 
                       private$fitting$fitted_result$coefficient[[1]]
                     compute_regularized_path_cpp(
                       fitting_train_i$reduced_data,
                       fitting_train_i$reduced_model,
                       fitting_train_i$control,
                       fitting_train_i$supplied_result,
                       fitting_train_i$fitted_result
                     )
                     subset_test_i <-
                       private$fitting$control$subset[private$fitting$control$cv_idx == i]
                     control$subset <- subset_test_i
                     fitting_test_i <-
                       lslxFitting$new(model = private$model,
                                       data = private$data,
                                       control = control)
                     cv_error_i <-
                       sapply(
                         X = 1:length(name_grid),
                         FUN = function(j) {
                           cv_error_ij <-
                             compute_loss_value_cpp(
                               theta_value = fitting_train_i$fitted_result$coefficient[[j]],
                               reduced_data = fitting_test_i$reduced_data,
                               reduced_model = fitting_test_i$reduced_model,
                               control = fitting_test_i$control,
                               supplied_result = fitting_test_i$supplied_result)
                           return(cv_error_ij)
                         }
                       )
                   }
                 )
               private$fitting$fitted_result$cv_error <- 
                 lapply(X = apply(cv_error, 1, mean),
                        FUN = function(cv_error_i) {
                          names(cv_error_i) <- "test_loss"
                          return(cv_error_i)
                        })
             }
           }
           
           names(private$fitting$fitted_result$numerical_condition) <-
             name_grid
           names(private$fitting$fitted_result$information_criterion) <-
             name_grid
           names(private$fitting$fitted_result$fit_index) <-
             name_grid
           names(private$fitting$fitted_result$cv_error) <-
             name_grid
           names(private$fitting$fitted_result$coefficient) <-
             name_grid
           
           idc_problem <-
             sapply(
               X = private$fitting$fitted_result$numerical_condition,
               FUN = function(numerical_condition_i) {
                 idc_problem_i <-
                   (numerical_condition_i[["n_iter_out"]] == private$fitting$control$iter_out_max) &
                   (numerical_condition_i[["objective_gradient_abs_max"]] > private$fitting$control$tol_out)
                 return(idc_problem_i)
               }
             )
           
           if (all(is.na(idc_problem))) {
             stop("Optimization result is incorrect under all specified penalty levels.\n",
                  "  Please check model identifiability or try other optimization parameters.\n")
           }
           
           if (verbose) {
             if (any(idc_problem)) {
               cat("WARNING: The algorithm may not converge under some penalty level. ")
               cat("Please try larger value of 'iter_out_max' or specify better starting values. \n")
             } else {
               cat("CONGRATS: The algorithm converged under all specified penalty levels. \n")
               cat("  Specified Tolerance for Convergence:",
                   private$fitting$control$tol_out,
                   "\n")
               cat(
                 "  Specified Maximal Number of Iterations:",
                 private$fitting$control$iter_out_max,
                 "\n"
               )
             }
           }
         })

## \code{$fit_lasso()} fits the specified model to data by minimizing a ML loss function with lasso penalty (Tibshirani, 1996). ##
lslx$set("public",
         "fit_lasso",
         function(lambda_grid = 0,
                  ...) {
           self$fit(penalty_method = "lasso",
                    lambda_grid = lambda_grid,
                    ...)
         })

## \code{$fit_mcp()} method fits the specified model to data by minimizing a ML loss function with mcp (Zhang, 2010). ##
lslx$set("public",
         "fit_mcp",
         function(lambda_grid = 0,
                  delta_grid = "default",
                  ...) {
           self$fit(
             penalty_method = "mcp",
             lambda_grid = lambda_grid,
             delta_grid = delta_grid,
             ...
           )
         })