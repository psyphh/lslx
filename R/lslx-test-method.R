lslx$set("public",
         "test_likelihood_ratio",
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
           
           fitted_numerical_condition <-
             self$extract_numerical_condition(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )[[1]]
           
           fitted_goodness_of_fit <-
             self$extract_goodness_of_fit(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )[[1]]
           
           likelihood_ratio <- list()
           likelihood_ratio$statistic <- 
             fitted_goodness_of_fit[["loss"]] * private$fitting$reduced_data$total_sample_size
           likelihood_ratio$degree_of_freedom <- 
             fitted_numerical_condition[["degree_of_freedom"]]
           likelihood_ratio$p_value <-
             1 - pchisq(likelihood_ratio$statistic,
                        likelihood_ratio$degree_of_freedom)
           return(likelihood_ratio)
         })


lslx$set("public",
         "test_rmsea",
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
           
           fitted_numerical_condition <-
             self$extract_numerical_condition(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )[[1]]
           
           fitted_goodness_of_fit <-
             self$extract_goodness_of_fit(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )[[1]]
           
           likelihood_ratio_statistic <-
             fitted_goodness_of_fit[["loss"]] * private$fitting$reduced_data$total_sample_size
           
           degree_of_freedom <-
             fitted_numerical_condition[["degree_of_freedom"]]
           
           rmsea <- list()
           rmsea$rmsea <- fitted_goodness_of_fit[["rmsea"]]
           
           lower_ncp <- 0
           if (pchisq(likelihood_ratio_statistic, degree_of_freedom, lower_ncp) < (1 - .05 / 2)) {
             
           } else {
             lower_ncp_1 <- lower_ncp
             lower_ncp_2 <- 0
             while (pchisq(likelihood_ratio_statistic,
                           degree_of_freedom,
                           lower_ncp_2) > (1 - .05 / 2)) {
               lower_ncp_2 <- lower_ncp_2 + degree_of_freedom
             }
             lower_ncp <- (lower_ncp_1 + lower_ncp_2) / 2
             while (abs(pchisq(likelihood_ratio_statistic, degree_of_freedom, lower_ncp) - (1 - .05 / 2)) > 10e-7) {
               if (pchisq(likelihood_ratio_statistic,
                          degree_of_freedom,
                          lower_ncp) < (1 - .05 / 2)) {
                 lower_ncp_2 <- lower_ncp
                 lower_ncp <- (lower_ncp + lower_ncp_1) / 2
               } else {
                 lower_ncp_1 <- lower_ncp
                 lower_ncp <- (lower_ncp + lower_ncp_2) / 2
                 
               }
             }
           }
           upper_ncp <- 0
           if (pchisq(likelihood_ratio_statistic, degree_of_freedom, upper_ncp) < (.05 / 2)) {
             upper_ncp <- 0
           } else {
             upper_ncp_1 <- upper_ncp
             upper_ncp_2 <- 0
             while (pchisq(likelihood_ratio_statistic,
                           degree_of_freedom,
                           upper_ncp_2) > (.05 / 2)) {
               upper_ncp_2 <- upper_ncp_2 + degree_of_freedom
             }
             upper_ncp <- (upper_ncp_1 + upper_ncp_2) / 2
             while (abs(pchisq(likelihood_ratio_statistic, degree_of_freedom, upper_ncp) - (.05 / 2)) > 10e-7) {
               if (pchisq(likelihood_ratio_statistic,
                          degree_of_freedom,
                          upper_ncp) < (.05 / 2)) {
                 upper_ncp_2 <- upper_ncp
                 upper_ncp <- (upper_ncp + upper_ncp_1) / 2
               } else {
                 upper_ncp_1 <- upper_ncp
                 upper_ncp <- (upper_ncp + upper_ncp_2) / 2
               }
             }
           }
           
           rmsea$lower_limit <-
             sqrt(
               private$fitting$reduced_model$n_group * lower_ncp /
                 (
                   private$fitting$reduced_data$total_sample_size * degree_of_freedom
                 )
             )
           rmsea$upper_limit <-
             sqrt(
               private$fitting$reduced_model$n_group * upper_ncp /
                 (
                   private$fitting$reduced_data$total_sample_size * degree_of_freedom
                 )
             )
           return(rmsea)
         })



lslx$set("public",
         "test_coefficient",
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
           
           fitted_coefficient <-
             self$extract_coefficient(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )[[1]]
           
           coefficient_acov <-
             self$extract_coefficient_acov(
               selector = selector,
               standard_error = standard_error,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           
           standard_error <- sqrt(diag(coefficient_acov))
           
           coefficient <- list()
           coefficient$estimate <- fitted_coefficient
           coefficient$standard_error <- standard_error
           coefficient$z_value <- fitted_coefficient / standard_error
           coefficient$p_value <- pnorm(-abs(coefficient$z_value))
           return(coefficient)
         })
