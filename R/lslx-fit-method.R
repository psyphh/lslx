lslx$set("public",
         "fit_mcp",
         function(lambda_grid = 0,
                  gamma_grid = Inf,
                  ...) {
           self$fit(
             penalty_method = "mcp",
             lambda_grid = lambda_grid,
             gamma_grid = gamma_grid,
             ...
           )
         })

lslx$set("public",
         "fit_lasso",
         function(lambda_grid = 0,
                  ...) {
           self$fit(penalty_method = "lasso",
                    lambda_grid = lambda_grid,
                    ...)
         })

lslx$set("public",
         "fit",
         function(penalty_method = "none",
                  lambda_grid = 0,
                  gamma_grid = Inf,
                  start_method = "MH",
                  positive_diag = TRUE,
                  iter_out_max = 100L,
                  iter_in_max = 50L,
                  iter_armijo_max = 100L,
                  tol_out = 1e-4,
                  tol_in = 1e-4,
                  step_size = 0.5,
                  armijo = 1e-5,
                  ridge_cov = 1e-4,
                  ridge_hessian = 1e-4,
                  verbose = TRUE) {
           if (!(penalty_method %in% c("none", "lasso", "mcp"))) {
             stop("Argument 'penalty_method' can be only either 'none', 'lasso', or 'mcp'.")
           } else {
             if (penalty_method == "none") {
               lambda_grid <- 0
               gamma_grid <- Inf
             }
             if (penalty_method == "lasso") {
               if (any(lambda_grid < 0)) {
                 stop(
                   "When argument 'penalty_method' is set as 'lasso', any element in argument 'lambda_grid' cannot be smaller than 0."
                 )
               }
               gamma_grid <- Inf
             }
             if (penalty_method == "mcp") {
               if (any(lambda_grid < 0)) {
                 stop(
                   "When argument 'penalty_method' is set as 'mcp', any element in argument 'lambda_grid' cannot be smaller than 0."
                 )
               }
               if (any(gamma_grid < 1)) {
                 stop(
                   "When argument 'penalty_method' is set as 'mcp', any element in argument 'gamma_grid' cannot be smaller than 1."
                 )
               }
             }
           }
           if (!(start_method %in% c("MH", "heuristic"))) {
             stop("Argument 'start_method' can be only 'MH' or 'heuristic'")
           }
           if (!(is.logical(positive_diag) &
                 (length(positive_diag) = 1))) {
             stop("Argument 'positive_diag' must be a numeric vector with length one.")
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
           
           control <-
             list(
               penalty_method = penalty_method,
               lambda_grid = lambda_grid,
               gamma_grid = gamma_grid,
               start_method = start_method,
               positive_diag = positive_diag,
               iter_out_max = iter_out_max,
               iter_in_max = iter_in_max,
               iter_armijo_max = iter_armijo_max,
               tol_out = tol_out,
               tol_in = tol_in,
               step_size = step_size,
               armijo = armijo,
               ridge_cov = ridge_cov,
               ridge_hessian = ridge_hessian
             )
           
           private$fitting <-
             lslxFitting$new(model = private$model,
                             data = private$data,
                             control = control)

           compute_regularized_path_cpp(
             private$fitting$reduced_data,
             private$fitting$reduced_model,
             private$fitting$control,
             private$fitting$numerical_condition,
             private$fitting$goodness_of_fit,
             private$fitting$coefficient
           )
           
           name_grid <-
             paste0(
               "ld=",
               sapply(
                 X = private$fitting$numerical_condition,
                 FUN = function(x) {
                   getElement(x, "lambda")
                 }
               ),
               "/",
               "gm=",
               sapply(
                 X = private$fitting$numerical_condition,
                 FUN = function(x) {
                   getElement(x, "gamma")
                 }
               )
             )
           
           names(private$fitting$numerical_condition) <- name_grid
           names(private$fitting$goodness_of_fit) <- name_grid
           names(private$fitting$coefficient) <- name_grid
         })
