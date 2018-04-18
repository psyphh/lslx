## \code{$validate()} prints a summary for the validation result under the given selector.. ##
lslx$set("public",
         "validate",
         function(data,
                  selector,
                  do_fit = "default",
                  standard_error = "default",
                  alpha_level = .05,
                  exclude_improper = TRUE) {
           if (missing(data)) {
             stop("Argument 'data' cannot be empty.")
           }
           if (is.null(private$fitting)) {
             stop("Fitting field is not yet derived. Please use fit-related methods first.")
           }
           if (!(
             do_fit %in% c("default", "none", "level", "pattern")
           )) {
             stop(
               "Argument 'do_fit' can be only either 'default', 'none', 'level', or 'pattern'."
             )
           }
           if (!(
             standard_error %in% c("default", "sandwich", "observed_fisher", "expected_fisher")
           )) {
             stop(
               "Argument 'standard_error' can be only either 'default', 'sandwich', 'observed_fisher', or 'expected_fisher'."
             )
           }
           if (do_fit == "default") {
             do_fit <- "pattern"
           }
           if (standard_error == "default") {
             if (private$fitting$control$response) {
               standard_error <- "sandwich"
             } else {
               standard_error <- "observed_fisher"
             }
           }
           
           lslx_cv <- self$clone(deep = TRUE)
           lslx_cv$set_data(data = data,
                            sample_cov = sample_cov,
                            sample_mean = sample_mean,
                            sample_size = sample_size)
           if (do_fit == "none") {
             
           } else {
             coefficient <-
               self$extract_coefficient(selector = selector,
                                        exclude_improper = exclude_improper)
             if (do_fit == "level") {
               penalty_level <-
                 self$extract_penalty_level(selector = selector,
                                            exclude_improper = exclude_improper)
               lambda <- as.numeric(strsplit(x = penalty_level,
                                             split = "=|/")[[1]][2])
               delta <- as.numeric(strsplit(x = penalty_level,
                                            split = "=|/")[[1]][4])
               lslx_cv$set_coefficient_start(name = names(coefficient), 
                                             start = coefficient)
               lslx_cv$fit(penalty_method = private$fitting$control$penalty_method,
                           lambda_grid = lambda,
                           delta_grid = delta,
                           algorithm = private$fitting$control$algorithm,
                           missing_method = private$fitting$control$missing_method,
                           start_method = "none",
                           iter_out_max = private$fitting$control$iter_out_max,
                           iter_in_max = private$fitting$control$iter_in_max,
                           iter_other_max = private$fitting$control$iter_other_max,
                           iter_armijo_max = private$fitting$control$iter_armijo_max,
                           tol_out = private$fitting$control$tol_out,
                           tol_in = private$fitting$control$tol_in,
                           tol_other = private$fitting$control$tol_other,
                           step_size = private$fitting$control$step_size,
                           armijo = private$fitting$control$armijo,
                           ridge_cov = private$fitting$control$ridge_cov,
                           ridge_hessian = private$fitting$control$ridge_hessian,
                           positive_diag = private$fitting$control$positive_diag,
                           verbose = FALSE)
             } else if (do_fit == "pattern") {
               type <- ifelse((coefficient != 0) & (private$model$specification$type != "fixed"),
                              "free", "fixed")
               lslx_cv$set_coefficient_start(name = names(coefficient), 
                                             start = coefficient)
               lslx_cv$set_coefficient_type(name = names(coefficient), 
                                            type = type)
               lslx_cv$fit(penalty_method = "none",
                           algorithm = private$fitting$control$algorithm,
                           missing_method = private$fitting$control$missing_method,
                           start_method = "none",
                           iter_out_max = private$fitting$control$iter_out_max,
                           iter_in_max = private$fitting$control$iter_in_max,
                           iter_other_max = private$fitting$control$iter_other_max,
                           iter_armijo_max = private$fitting$control$iter_armijo_max,
                           tol_out = private$fitting$control$tol_out,
                           tol_in = private$fitting$control$tol_in,
                           tol_other = private$fitting$control$tol_other,
                           step_size = private$fitting$control$step_size,
                           armijo = private$fitting$control$armijo,
                           ridge_cov = private$fitting$control$ridge_cov,
                           ridge_hessian = private$fitting$control$ridge_hessian,
                           positive_diag = private$fitting$control$positive_diag,
                           verbose = FALSE)
             } else {
             }
             lslx_cv$summarize(selector = selector,
                               standard_error = standard_error,
                               alpha_level = alpha_level,
                               exclude_improper = exclude_improper)
           } 
         })
