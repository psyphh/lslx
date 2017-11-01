lslx$set("public",
         "extract_specification",
         function() {
           specification <-
             private$model$specification
           return(specification)
         })


lslx$set("public",
         "extract_saturated_cov",
         function() {
           saturated_cov <-
             private$fitting$reduced_data$saturated_cov
           return(saturated_cov)
         })


lslx$set("public",
         "extract_saturated_mean",
         function() {
           saturated_mean <-
             private$fitting$reduced_data$saturated_mean
           return(saturated_mean)
         })


lslx$set("public",
         "extract_saturated_moment_acov",
         function() {
           saturated_moment_acov <-
             private$fitting$reduced_data$saturated_moment_acov
           return(saturated_moment_acov)
         })


lslx$set("public",
         "extract_penalty_level",
         function(selector,
                  exclude_improper = TRUE) {
           if (missing(selector)) {
             if (length(private$fitting$fitted_result$numerical_condition) == 1) {
               selector <- "bic"
             } else {
               stop("Argument 'selector' is missing.")
             }
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           
           if (!(selector %in% names(private$fitting$fitted_result$information_criterion[[1]]))) {
             stop(
               "Argument 'selector' is unrecognized.",
               "\n  Selector currently recognized by 'lslx' is \n  ",
               do.call(paste, as.list(
                 names(private$fitting$fitted_result$information_criterion[[1]])
               )),
               "."
             )
           }
           if ((selector %in% c("raic", "raic3", "rcaic", "rbic", "rabic", "rhbic")) & 
               (!private$fitting$control$response)) {
             stop("When lslx object is initialized via moments,",
                  " 'raic', 'raic3', 'rcaic', 'rbic', 'rabic', and 'rhbic' are not available.")
           }
           
           if (exclude_improper) {
             idx_convergence <-
               which(
                 sapply(
                   X = private$fitting$fitted_result$numerical_condition,
                   FUN = function(x) {
                     getElement(object = x,
                                name = "objective_gradient_abs_max")
                   }
                 ) < private$fitting$control$tol_out
               )
             if (length(idx_convergence) == 0) {
               stop(
                 "The PL estimates under all penalty levels are derived under nonconverged result. \n",
                 "Please try larger value of 'iter_out_max' or specify better starting values."
               )
             }
             idx_convexity <-
               which(
                 sapply(
                   X = private$fitting$fitted_result$numerical_condition,
                   FUN = function(x) {
                     getElement(object = x,
                                name = "objective_hessian_convexity")
                   }
                 ) > 0
               )
             if (length(idx_convexity) == 0) {
               stop(
                 "The PL estimates under all penalty levels are derived under nonconvex hessian. \n",
                 "Please try larger value of delta or check the identifiability of specified model."
               )
             }
           } else {
             idx_convergence <-
               seq_len(length(private$fitting$fitted_result$numerical_condition))
             idx_convexity <-
               seq_len(length(private$fitting$fitted_result$numerical_condition))
           }
           idx_selection <-
             intersect(x = idx_convergence, y = idx_convexity)
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
                   names(which.min(information_criterion_i[idx_selection]))
                 return(penalty_level_i)
               }
             )
           return(penalty_level)
         })


lslx$set("public",
         "extract_numerical_condition",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           numerical_condition <-
             private$fitting$fitted_result$numerical_condition[[penalty_level]]
           return(numerical_condition)
         })


lslx$set("public",
         "extract_information_criterion",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           information_criterion <-
             private$fitting$fitted_result$information_criterion[[penalty_level]]
           return(information_criterion)
         })

lslx$set("public",
         "extract_fit_indice",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           fit_indice <-
             private$fitting$fitted_result$fit_indice[[penalty_level]]
           return(fit_indice)
         })


lslx$set("public",
         "extract_coefficient",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           return(coefficient)
         })


lslx$set("public",
         "extract_implied_cov",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
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
           names(implied_cov) <- private$model$name_group
           return(implied_cov)
         })

lslx$set("public",
         "extract_implied_mean",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
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
           names(implied_mean) <- private$model$name_group
           return(implied_mean)
         })


lslx$set("public",
         "extract_residual_cov",
         function(selector,
                  exclude_improper = TRUE) {
           implied_cov <-
             self$extract_implied_cov(selector = selector,
                                      exclude_improper = exclude_improper)
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

lslx$set("public",
         "extract_residual_mean",
         function(selector,
                  exclude_improper = TRUE) {
           implied_mean <-
             self$extract_implied_mean(selector = selector,
                                       exclude_improper = exclude_improper)
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

lslx$set("public",
         "extract_coefficient_matrice",
         function(selector,
                  block,
                  exclude_improper = TRUE) {
           if (missing(block)) {
             stop("Argument 'block' is missing.")
           }
           if (length(block) > 1) {
             stop("The length of argument 'block' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           coefficient_matrice <-
             compute_coefficient_matrice_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           
           if (block %in% c("f<-1", "y<-1")) {
             selected_matrix <- coefficient_matrice$alpha
             col_names <- "1"
           } else {
             if (block %in% c("f<-f", "f<-y", "y<-f", "y<-y")) {
               selected_matrix <- coefficient_matrice$beta
             } else if (block %in% c("f<->f", "f<->y", "y<->f", "y<->y")) {
               selected_matrix <- coefficient_matrice$psi
             } else {
               stop(
                 "Argument 'block' is unrecognized. It must be one of the following:\n",
                 "  'f<-1', 'y<-1', 'f<-f', 'f<-y', 'y<-f', 'f<-f', 'f<->f', 'f<->y', 'y<->f', 'f<->f'."
               )
             }
             col_names <- private$model$name_eta
           }
           row_select <- strsplit(block, split = "<-|<->|->")[[1]][1]
           col_select <- strsplit(block, split = "<-|<->|->")[[1]][2]
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
           coefficient_matrice_block <-
             lapply(
               X = selected_matrix,
               FUN = function(selected_matrix_i) {
                 coefficient_matrice_block_i <-
                   selected_matrix_i[row_select,
                                     col_select,
                                     drop = FALSE]
                 return(coefficient_matrice_block_i)
               }
             )
           names(coefficient_matrice_block) <-
             private$model$name_group
           return(coefficient_matrice_block)
         })




lslx$set("public",
         "extract_moment_jacobian",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           moment_jacobian <-
             compute_moment_jacobian_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(moment_jacobian) <-
             rownames(private$model$specification)
           return(moment_jacobian)
         })



lslx$set("public",
         "extract_expected_fisher",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           expected_fisher <-
             compute_expected_fisher_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(expected_fisher) <-
             rownames(private$model$specification)
           rownames(expected_fisher) <-
             rownames(private$model$specification)
           return(expected_fisher)
         })


lslx$set("public",
         "extract_observed_fisher",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           observed_fisher <-
             compute_observed_fisher_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           colnames(observed_fisher) <-
             rownames(private$model$specification)
           rownames(observed_fisher) <-
             rownames(private$model$specification)
           return(observed_fisher)
         })


lslx$set("public",
         "extract_score_acov",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
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
           return(score_acov)
         })


lslx$set("public",
         "extract_coefficient_acov",
         function(selector,
                  standard_error = "default",
                  exclude_improper = TRUE) {
           if (!(
             standard_error %in% c("default", "sandwich", "observed_fisher", "expected_fisher")
           )) {
             stop(
               "Argument 'standard_error' can be only either 'default', 'sandwich', 'observed_fisher', or 'expected_fisher'."
             )
           }
           if (standard_error == "default") {
             if (private$fitting$control$response) {
               standard_error <- "sandwich"
             } else {
               standard_error <- "observed_fisher"
             }
           }
           coefficient <-
             self$extract_coefficient(selector = selector,
                                      exclude_improper = exclude_improper)
           idc_effective <-
             private$fitting$reduced_model$theta_is_free |
             (private$fitting$reduced_model$theta_is_pen &
                coefficient != 0)
           coefficient_acov <-
             matrix(NA, length(coefficient), length(coefficient))
           colnames(coefficient_acov) <-
             rownames(private$model$specification)
           rownames(coefficient_acov) <-
             rownames(private$model$specification)
           if (standard_error == "sandwich") {
             score_acov <-
               self$extract_score_acov(selector = selector,
                                       exclude_improper = exclude_improper)
             observed_fisher <-
               self$extract_observed_fisher(selector = selector,
                                            exclude_improper = exclude_improper)
             observed_fisher_pinv <-
               solve(observed_fisher[idc_effective, idc_effective])
             coefficient_acov[idc_effective, idc_effective] <-
               (observed_fisher_pinv %*%
                  score_acov[idc_effective, idc_effective] %*%
                  observed_fisher_pinv)
           } else if (standard_error == "expected_fisher") {
             expected_fisher <-
               self$extract_expected_fisher(selector = selector,
                                            exclude_improper = exclude_improper)
             coefficient_acov[idc_effective, idc_effective] <-
               solve(expected_fisher[idc_effective, idc_effective]) /
               private$fitting$reduced_data$n_observation
           } else if (standard_error == "observed_fisher") {
             observed_fisher <-
               self$extract_observed_fisher(selector = selector,
                                            exclude_improper = exclude_improper)
             coefficient_acov[idc_effective, idc_effective] <-
               solve(observed_fisher[idc_effective, idc_effective]) /
               private$fitting$reduced_data$n_observation
           } else {
           }
           attr(coefficient_acov, "standard_error") <- standard_error
           return(coefficient_acov)
         })




lslx$set("public",
         "extract_loss_gradient",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           loss_gradient <-
             compute_loss_gradient_direct_cpp(
               theta_value = coefficient,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           rownames(loss_gradient) <-
             rownames(private$model$specification)
           return(loss_gradient)
         })



lslx$set("public",
         "extract_regularizer_gradient",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           lambda <- as.numeric(strsplit(x = penalty_level,
                                         split = "=|/")[[1]][2])
           delta <- as.numeric(strsplit(x = penalty_level,
                                        split = "=|/")[[1]][4])
           regularizer_gradient <-
             compute_regularizer_gradient_cpp(
               theta_value = coefficient,
               lambda = lambda,
               delta = delta,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           rownames(regularizer_gradient) <-
             rownames(private$model$specification)
           return(regularizer_gradient)
         })

lslx$set("public",
         "extract_objective_gradient",
         function(selector,
                  exclude_improper = TRUE) {
           penalty_level <-
             self$extract_penalty_level(selector = selector,
                                        exclude_improper = exclude_improper)
           coefficient <-
             private$fitting$fitted_result$coefficient[[penalty_level]]
           lambda <- as.numeric(strsplit(x = penalty_level,
                                         split = "=|/")[[1]][2])
           delta <- as.numeric(strsplit(x = penalty_level,
                                        split = "=|/")[[1]][4])
           objective_gradient <-
             compute_objective_gradient_cpp(
               theta_value = coefficient,
               lambda = lambda,
               delta = delta,
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               supplied_result = private$fitting$supplied_result
             )
           rownames(objective_gradient) <-
             rownames(private$model$specification)
           return(objective_gradient)
         })
