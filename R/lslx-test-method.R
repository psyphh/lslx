## \code{$test_lr()} returns a \code{data.frame} of result for likelihood ratio test. ##
lslx$set("public",
         "test_lr",
         function(selector,
                  exclude_improper = TRUE) {
           numerical_condition <-
             self$extract_numerical_condition(selector = selector,
                                              exclude_improper = exclude_improper)
           lr_test <-
             data.frame(
               statistic = c(NA, NA),
               df = c(NA, NA),
               p_value = c(NA, NA),
               row.names = c("unadjusted", "mean-adjusted")
             )
           lr_test["unadjusted", "statistic"] <-
             numerical_condition[["loss_value"]] * private$fitting$reduced_data$n_observation
           lr_test["unadjusted", "df"] <-
             numerical_condition[["degree_of_freedom"]]
           lr_test["unadjusted", "p_value"] <-
             1 - pchisq(lr_test["unadjusted", "statistic"],
                        lr_test["unadjusted", "df"])
           if (private$fitting$control$response) {
             scaling_factor <- numerical_condition[["scaling_factor"]]
             if (!is.na(scaling_factor)) {
               lr_test["mean-adjusted", "statistic"] <-
                 numerical_condition[["loss_value"]] * private$fitting$reduced_data$n_observation /
                 scaling_factor
               lr_test["mean-adjusted", "df"] <-
                 numerical_condition[["degree_of_freedom"]]
               lr_test["mean-adjusted", "p_value"] <-
                 1 - pchisq(lr_test["mean-adjusted", "statistic"],
                            lr_test["mean-adjusted", "df"])
             }
           } else {
             
           }
           return(lr_test)
         })

## \code{$test_rmsea()} returns a \code{data.frame} of result for rmsea confidence intervals. ##
lslx$set("public",
         "test_rmsea",
         function(selector,
                  alpha_level = .05,
                  exclude_improper = TRUE) {
           numerical_condition <-
             self$extract_numerical_condition(selector = selector,
                                              exclude_improper = exclude_improper)
           fit_indice <-
             self$extract_fit_indice(selector = selector,
                                     exclude_improper = exclude_improper)
           lr_test <-
             self$test_lr(selector = selector,
                          exclude_improper = exclude_improper)
           rmsea_test <-
             data.frame(
               estimate = c(NA, NA),
               lower = c(NA, NA),
               upper = c(NA, NA),
               row.names = c("unadjusted", "mean-adjusted")
             )
           for (row_name_i in row.names(rmsea_test)) {
             if ((row_name_i == "unadjusted") |
                 (private$fitting$control$response)) {
               lr_statistic <-
                 lr_test[row_name_i, "statistic"]
               lr_df <- lr_test[row_name_i, "df"]
               if (is.na(lr_statistic) | is.na(lr_df)) {
                 rmsea_test[row_name_i, "estimate"] <- NA
                 rmsea_test[row_name_i, "lower"] <- NA
                 rmsea_test[row_name_i, "upper"] <- NA
               } else if ((lr_df == 0) & (lr_statistic > sqrt(.Machine$double.eps)) ) {
                 rmsea_test[row_name_i, "estimate"] <- NA
                 rmsea_test[row_name_i, "lower"] <- NA
                 rmsea_test[row_name_i, "upper"] <- NA
               } else if (lr_statistic < sqrt(.Machine$double.eps)) {
                 rmsea_test[row_name_i, "estimate"] <- 0
                 rmsea_test[row_name_i, "lower"] <- 0
                 rmsea_test[row_name_i, "upper"] <- 0
               } else {
                 lower_ncp <- 0
                 if (pchisq(lr_statistic,
                            lr_df, lower_ncp) < (1 - alpha_level / 2)) {
                 } else {
                   lower_ncp_1 <- lower_ncp
                   lower_ncp_2 <- 0
                   while (pchisq(lr_statistic,
                                 lr_df,
                                 lower_ncp_2) > (1 - alpha_level / 2)) {
                     lower_ncp_2 <- lower_ncp_2 + lr_df
                   }
                   lower_ncp <- (lower_ncp_1 + lower_ncp_2) / 2
                   while (abs(pchisq(lr_statistic,
                                     lr_df, lower_ncp) -
                              (1 - alpha_level / 2)) > private$fitting$control$tol_other) {
                     if (pchisq(lr_statistic,
                                lr_df,
                                lower_ncp) < (1 - alpha_level / 2)) {
                       lower_ncp_2 <- lower_ncp
                       lower_ncp <- (lower_ncp + lower_ncp_1) / 2
                     } else {
                       lower_ncp_1 <- lower_ncp
                       lower_ncp <- (lower_ncp + lower_ncp_2) / 2
                     }
                   }
                 }
                 upper_ncp <- 0
                 if (pchisq(lr_statistic,
                            lr_df, upper_ncp) < (alpha_level / 2)) {
                 } else {
                   upper_ncp_1 <- upper_ncp
                   upper_ncp_2 <- 0
                   while (pchisq(lr_statistic,
                                 lr_df,
                                 upper_ncp_2) > (alpha_level / 2)) {
                     upper_ncp_2 <- upper_ncp_2 + lr_df
                   }
                   upper_ncp <- (upper_ncp_1 + upper_ncp_2) / 2
                   while (abs(pchisq(lr_statistic,
                                     lr_df, upper_ncp) -
                              (alpha_level / 2)) > private$fitting$control$tol_other) {
                     if (pchisq(lr_statistic,
                                lr_df,
                                upper_ncp) < (alpha_level / 2)) {
                       upper_ncp_2 <- upper_ncp
                       upper_ncp <- (upper_ncp + upper_ncp_1) / 2
                     } else {
                       upper_ncp_1 <- upper_ncp
                       upper_ncp <- (upper_ncp + upper_ncp_2) / 2
                     }
                   }
                 }
                 if (row_name_i == "unadjusted") {
                   scaling_factor <- 1
                 } else {
                   scaling_factor <- numerical_condition[["scaling_factor"]]
                 }
                 if (!is.na(scaling_factor)) {
                   rmsea_test[row_name_i, "estimate"] <-
                     sqrt(max(
                       0,
                       scaling_factor * private$fitting$reduced_model$n_group * (lr_statistic - lr_df) /
                         (private$fitting$reduced_data$n_observation * lr_df)
                     ))
                   rmsea_test[row_name_i, "lower"]  <-
                     sqrt(max(
                       0,
                       scaling_factor * private$fitting$reduced_model$n_group * lower_ncp /
                         (private$fitting$reduced_data$n_observation * lr_df)
                     ))
                   rmsea_test[row_name_i, "upper"]  <-
                     sqrt(max(
                       0,
                       scaling_factor * private$fitting$reduced_model$n_group * upper_ncp /
                         (private$fitting$reduced_data$n_observation * lr_df)
                     ))
                 }
               }
             }
           }
           return(rmsea_test)
         })

## \code{$test_coefficient()} returns a \code{data.frame} of result for coefficient significance and confidence interval. ##
lslx$set("public",
         "test_coefficient",
         function(selector,
                  standard_error = "default",
                  alpha_level = .05,
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
           coefficient_acov <-
             self$extract_coefficient_acov(
               selector = selector,
               standard_error = standard_error,
               exclude_improper = exclude_improper
             )
           coefficient_test <-
             data.frame(estimate = coefficient,
                        standard_error = sqrt(diag(coefficient_acov)))
           coefficient_test$z_value <-
             coefficient_test$estimate / coefficient_test$standard_error
           coefficient_test$p_value <-
             pnorm(-abs(coefficient_test$z_value))
           coefficient_test$lower <-
             coefficient_test$estimate + qnorm(alpha_level / 2) * coefficient_test$standard_error
           coefficient_test$upper <-
             coefficient_test$estimate + qnorm(1 - alpha_level / 2) * coefficient_test$standard_error
           attr(coefficient_test, "standard_error") <-
             standard_error
           return(coefficient_test)
         })
