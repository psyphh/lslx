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
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (!all(selector %in% names(private$fitting$goodness_of_fit[[1]]))) {
             stop(
               "Argument 'selector' contains unrecognized fit indice(s).",
               "\n  Fit indice(s) currently recognized by 'lslx' is \n  ",
               do.call(paste, as.list(
                 names(private$fitting$goodness_of_fit[[1]])
               )),
               "."
             )
           }
           if (exclude_nonconvergence) {
             idx_convergence <-
               which(
                 sapply(
                   X = private$fitting$numerical_condition,
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
           } else {
             idx_convergence <-
               seq_len(length(private$fitting$numerical_condition))
           }
           
           if (exclude_nonconvexity) {
             idx_convexity <-
               which(sapply(
                 X = private$fitting$numerical_condition,
                 FUN = function(x) {
                   getElement(object = x,
                              name = "objective_hessian_convexity")
                 }
               ) > 0)
             if (length(idx_convexity) == 0) {
               stop(
                 "The PL estimates under all penalty levels are derived under nonconvex hessian. \n",
                 "Please try larger value of gamma or check the identifiability of specified model."
               )
             }
           } else {
             idx_convexity <-
               seq_len(length(private$fitting$numerical_condition))
           }
           
           idx_selection <-
             intersect(x = idx_convergence, y = idx_convexity)
           
           penalty_level <-
             sapply(
               X = selector,
               FUN = function(selector_i) {
                 goodness_of_fit_i <- sapply(
                   X = private$fitting$goodness_of_fit,
                   FUN = function(goodness_of_fit_j) {
                     getElement(object = goodness_of_fit_j,
                                name = selector_i)
                   }
                 )
                 if (selector_i %in% c("cfi",
                                       "nnfi")) {
                   penalty_level_i <-
                     names(which.max(goodness_of_fit_i[idx_selection]))
                 } else {
                   penalty_level_i <-
                     names(which.min(goodness_of_fit_i[idx_selection]))
                 }
                 return(penalty_level_i)
               }
             )
           return(penalty_level)
         })



lslx$set("public",
         "extract_coefficient",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           coefficient <-
             private$fitting$coefficient[penalty_level]
           return(coefficient)
         })

lslx$set("public",
         "extract_goodness_of_fit",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           goodness_of_fit <-
             private$fitting$goodness_of_fit[penalty_level]
           return(goodness_of_fit)
         })


lslx$set("public",
         "extract_numerical_condition",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           numerical_condition <-
             private$fitting$numerical_condition[penalty_level]
           return(numerical_condition)
         })


lslx$set("public",
         "extract_implied_cov",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           implied_cov <-
             compute_implied_cov_cpp(
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               theta_value = fitted_coefficient
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
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           implied_mean <-
             compute_implied_mean_cpp(
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               theta_value = fitted_coefficient
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
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           implied_cov <-
             self$extract_implied_cov(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
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

lslx$set("public",
         "extract_residual_mean",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           implied_mean <-
             self$extract_implied_mean(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
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

lslx$set("public",
         "extract_coefficient_matrice",
         function(selector,
                  block,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           if (missing(block)) {
             stop("Argument 'block' is missing.")
           }
           if (length(block) > 1) {
             stop("The length of argument 'block' can be only one.")
           }
           
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           coefficient_matrice <-
             compute_coefficient_matrice_cpp(
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               theta_value = fitted_coefficient
             )
           if (block %in% c("f<-1", "y<-1")) {
             alpha <- coefficient_matrice$alpha
             alpha <-
               lapply(
                 X = alpha,
                 FUN = function(alpha_i) {
                   rownames(alpha_i) <- private$model$name_eta
                   colnames(alpha_i) <- "1"
                   return(alpha_i)
                 }
               )
             if (block == "f<-1") {
               coefficient_matrice_block <-
                 lapply(
                   X = alpha,
                   FUN = function(alpha_i) {
                     coefficient_matrice_block_i <-
                       alpha_i[private$model$name_factor, "1", drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             } else {
               coefficient_matrice_block <-
                 lapply(
                   X = alpha,
                   FUN = function(alpha_i) {
                     coefficient_matrice_block_i <-
                       alpha_i[private$model$name_response, "1", drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             }
           } else if (block %in% c("f<-f", "f<-y", "y<-f", "y<-y")) {
             beta <- coefficient_matrice$beta
             beta <-
               lapply(
                 X = beta,
                 FUN = function(beta_i) {
                   rownames(beta_i) <- private$model$name_eta
                   colnames(beta_i) <- private$model$name_eta
                   return(beta_i)
                 }
               )
             if (block == "f<-f") {
               coefficient_matrice_block <-
                 lapply(
                   X = beta,
                   FUN = function(beta_i) {
                     coefficient_matrice_block_i <-
                       beta_i[private$model$name_factor,
                              private$model$name_factor,
                              drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             } else if (block == "f<-y") {
               coefficient_matrice_block <-
                 lapply(
                   X = beta,
                   FUN = function(beta_i) {
                     coefficient_matrice_block_i <-
                       beta_i[private$model$name_factor,
                              private$model$name_response,
                              drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             } else if (block == "y<-f") {
               coefficient_matrice_block <-
                 lapply(
                   X = beta,
                   FUN = function(beta_i) {
                     coefficient_matrice_block_i <-
                       beta_i[private$model$name_response,
                              private$model$name_factor,
                              drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             } else {
               coefficient_matrice_block <-
                 lapply(
                   X = beta,
                   FUN = function(beta_i) {
                     coefficient_matrice_block_i <-
                       beta_i[private$model$name_response,
                              private$model$name_response,
                              drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             }
           } else if (block %in% c("f<->f", "f<->y", "y<->f", "y<->y")) {
             psi <- coefficient_matrice$psi
             psi <-
               lapply(
                 X = psi,
                 FUN = function(psi_i) {
                   rownames(psi_i) <- private$model$name_eta
                   colnames(psi_i) <- private$model$name_eta
                   return(psi_i)
                 }
               )
             if (block == "f<->f") {
               coefficient_matrice_block <-
                 lapply(
                   X = psi,
                   FUN = function(psi_i) {
                     coefficient_matrice_block_i <-
                       psi_i[private$model$name_factor,
                             private$model$name_factor,
                             drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             } else if (block == "f<->y") {
               coefficient_matrice_block <-
                 lapply(
                   X = psi,
                   FUN = function(psi_i) {
                     coefficient_matrice_block_i <-
                       psi_i[private$model$name_factor,
                             private$model$name_response,
                             drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             } else if (block == "y<->f") {
               coefficient_matrice_block <-
                 lapply(
                   X = psi,
                   FUN = function(psi_i) {
                     coefficient_matrice_block_i <-
                       psi_i[private$model$name_response,
                             private$model$name_factor,
                             drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             } else {
               coefficient_matrice_block <-
                 lapply(
                   X = psi,
                   FUN = function(psi_i) {
                     coefficient_matrice_block_i <-
                       psi_i[private$model$name_response,
                             private$model$name_response,
                             drop = FALSE]
                     return(coefficient_matrice_block_i)
                   }
                 )
             }
           } else {
             stop(
               "Argument 'block' is unrecognized. It must be one of the following:\n",
               "  'f<-1', 'y<-1', 'f<-f', 'f<-y', 'y<-f', 'f<-f', 'f<->f', 'f<->y', 'y<->f', 'f<->f'."
             )
           }
           names(coefficient_matrice_block) <-
             private$model$name_group
           return(coefficient_matrice_block)
         })




lslx$set("public",
         "extract_moment_gradient",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           moment_gradient <-
             compute_moment_gradient_cpp(
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               theta_value = fitted_coefficient
             )
           colnames(moment_gradient) <-
             rownames(private$model$specification)
           return(moment_gradient)
         })



lslx$set("public",
         "extract_fisher_information",
         function(selector,
                  type = "expected",
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           if (type == "expected") {
             fisher_information <-
               compute_expected_fisher_information_cpp(
                 reduced_data = private$fitting$reduced_data,
                 reduced_model = private$fitting$reduced_model,
                 control = private$fitting$control,
                 theta_value = fitted_coefficient
               )
           } else {
             stop("So far only expected Fisher information can be extracted.")
           }
           colnames(fisher_information) <-
             rownames(private$model$specification)
           rownames(fisher_information) <-
             rownames(private$model$specification)
           return(fisher_information)
         })

lslx$set("public",
         "extract_coefficient_acov",
         function(selector,
                  standard_error = "expected_fisher",
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           if (standard_error == "expected_fisher") {
             fisher_information <-
               self$extract_fisher_information(
                 selector = selector,
                 type = "expected",
                 exclude_nonconvergence = exclude_nonconvergence,
                 exclude_nonconvexity = exclude_nonconvexity,
                 verbose = verbose
               )
           } else {
             stop(
               "In the current version, only acov based on expected Fisher information can be extracted."
             )
           }
           fitted_coefficient <-
             self$extract_coefficient(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )[[1]]
           idc_effective <-
             private$fitting$reduced_model$theta_is_free |
             (private$fitting$reduced_model$theta_is_pen &
                fitted_coefficient != 0)
           coefficient_acov <-
             matrix(NA, length(fitted_coefficient), length(fitted_coefficient))
           colnames(coefficient_acov) <-
             rownames(private$model$specification)
           rownames(coefficient_acov) <-
             rownames(private$model$specification)
           
           coefficient_acov[idc_effective, idc_effective] <-
             solve(fisher_information[idc_effective, idc_effective]) /
             private$fitting$reduced_data$total_sample_size
           return(coefficient_acov)
         })





lslx$set("public",
         "extract_loss_gradient",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           loss_gradient <-
             compute_loss_gradient_cpp(
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               theta_value = fitted_coefficient
             )
           rownames(loss_gradient) <-
             rownames(private$model$specification)
           return(loss_gradient)
         })

lslx$set("public",
         "extract_regularizer_gradient",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           
           lambda <- as.numeric(strsplit(x = penalty_level,
                                         split = "=|/")[[1]][2])
           gamma <- as.numeric(strsplit(x = penalty_level,
                                        split = "=|/")[[1]][4])
           regularizer_gradient <-
             compute_regularizer_gradient_cpp(
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               theta_value = fitted_coefficient,
               lambda = lambda,
               gamma = gamma
             )
           rownames(regularizer_gradient) <-
             rownames(private$model$specification)
           return(regularizer_gradient)
         })

lslx$set("public",
         "extract_objective_gradient",
         function(selector,
                  exclude_nonconvergence = TRUE,
                  exclude_nonconvexity = TRUE,
                  verbose = TRUE) {
           if (missing(selector)) {
             stop("Argument 'selector' is missing.")
           }
           if (length(selector) > 1) {
             stop("The length of argument 'selector' can be only one.")
           }
           penalty_level <-
             self$extract_penalty_level(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           fitted_coefficient <-
             private$fitting$coefficient[[penalty_level]]
           lambda <- as.numeric(strsplit(x = penalty_level,
                                         split = "=|/")[[1]][2])
           gamma <- as.numeric(strsplit(x = penalty_level,
                                        split = "=|/")[[1]][4])
           objective_gradient <-
             compute_objective_gradient_cpp(
               reduced_data = private$fitting$reduced_data,
               reduced_model = private$fitting$reduced_model,
               control = private$fitting$control,
               theta_value = fitted_coefficient,
               lambda = lambda,
               gamma = gamma
             )
           rownames(objective_gradient) <-
             rownames(private$model$specification)
           return(objective_gradient)
         })
