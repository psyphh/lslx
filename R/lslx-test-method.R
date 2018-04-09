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
               } else if ((lr_df == 0) &
                          (lr_statistic > sqrt(.Machine$double.eps))) {
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
                     sqrt(
                       max(
                         0,
                         scaling_factor * private$fitting$reduced_model$n_group * lower_ncp /
                           (private$fitting$reduced_data$n_observation * lr_df)
                       )
                     )
                   rmsea_test[row_name_i, "upper"]  <-
                     sqrt(
                       max(
                         0,
                         scaling_factor * private$fitting$reduced_model$n_group * upper_ncp /
                           (private$fitting$reduced_data$n_observation * lr_df)
                       )
                     )
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
                  debias = "default",
                  post = "default",
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
           if (post == "default") {
             post <- "none"
             if (debias == "default") {
               debias <- FALSE
             }
           } else if (post == "polyhedral") {
             if (debias == "default") {
               debias <- TRUE
             }
             if (!debias) {
               stop(
                 "'debias' cannot be FALSE under 'post' == 'polyhedral'."
               )
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
           if (!debias) {
             coefficient_test <-
               data.frame(estimate = coefficient,
                          standard_error = sqrt(diag(coefficient_acov)))             
           } else {
             debiased_coefficient <-
               self$extract_debiased_coefficient(selector = selector,
                                                 exclude_improper = exclude_improper)
             coefficient_test <-
               data.frame(estimate = debiased_coefficient,
                          standard_error = sqrt(diag(coefficient_acov)))
           }
           coefficient_test$z_value <-
             coefficient_test$estimate / coefficient_test$standard_error
           if (post == "none") {
             coefficient_test$p_value <-
               pnorm(-abs(coefficient_test$z_value))
             coefficient_test$lower <-
               coefficient_test$estimate + qnorm(alpha_level / 2) * coefficient_test$standard_error
             coefficient_test$upper <-
               coefficient_test$estimate + qnorm(1 - alpha_level / 2) * coefficient_test$standard_error
             attr(coefficient_test, "standard_error") <-
               standard_error
           } else if (post == "polyhedral") {
             is_active <-
               private$fitting$reduced_model$theta_is_free |
               (private$fitting$reduced_model$theta_is_pen &
                  coefficient != 0)
             is_pen <- private$fitting$reduced_model$theta_is_pen
             is_selected <- is_pen & (coefficient != 0)
             if (!any(is_selected)) {
               stop(
                 "No non-zero parameters are selected and hence selective inference cannot be implemented."
               )
             }
             a_ph <- - diag(sign(coefficient))
             b_ph <-
               matrix((sign(coefficient) * (coefficient - debiased_coefficient)))
             coefficient_test$p_value <-
               sapply(
                 X = 1:length(debiased_coefficient),
                 FUN = function(i) {
                   tnorm_quantity <- 
                     compute_tnorm_quantity(i, a_ph, b_ph, 
                                         debiased_coefficient, coefficient_acov,
                                         is_pen, is_active, is_selected) 
                   tnorm_p_value <-  
                     compute_tnorm_p_value(theta = tnorm_quantity$theta, 
                                           mu = 0, sigma = tnorm_quantity$sigma, 
                                           left = tnorm_quantity$left, 
                                           right = tnorm_quantity$right)
                   return(tnorm_p_value)
                 },
                 simplify = TRUE,
                 USE.NAMES = TRUE
               )
             coefficient_test$lower <- NA_real_
             coefficient_test$upper <- NA_real_
           } else {
             
           }
           return(coefficient_test)
         })



compute_tnorm_quantity <- function(i, a_ph, b_ph, 
                                debiased_coefficient, coefficient_acov,
                                is_pen, is_active, is_selected) {
  theta_i <- debiased_coefficient[i]
  sigma_i <- sqrt(coefficient_acov[i, i])
  if (is_pen[i]) {
    if (is_selected[i]) {
      coefficient_acov[is.na(coefficient_acov)] <- 0
      a_ph <- a_ph[is_selected, is_active, drop = FALSE]
      b_ph <- b_ph[is_selected, 1, drop = FALSE]
      n_is_active <- sum(is_active)
      cumsum_is_active <- cumsum(is_active)
      c_i <- (coefficient_acov[is_active, i, drop = FALSE]) / (sigma_i^2)
      c_i_expand <- diag(0, n_is_active)
      c_i_expand[, cumsum_is_active[i]] <- c_i
      z_i <-
        (diag(n_is_active) - c_i_expand) %*% matrix(debiased_coefficient[is_active])
      w_i <- c(a_ph %*% c_i)
      r_i <- c(b_ph - (a_ph %*% z_i)) / w_i
      left_i <- ifelse(any(w_i < 0), yes = max(r_i[w_i < 0]), no = - Inf) 
      right_i <- ifelse(any(w_i > 0), yes = min(r_i[w_i > 0]), no = Inf) 
    } else {
      left_i <- NA
      right_i <- NA
    }
  } else {
    left_i <- - Inf
    right_i <- Inf
  }
  tnorm_quantity <- 
    list(theta = theta_i, sigma = sigma_i, left = left_i, right = right_i)
  return(tnorm_quantity)
}

compute_tnorm_p_value <- function(theta, mu, sigma, left, right) {
  if (is.na(left) | is.na(right) | is.na(sigma)) {
    p_value <- NA
  } else {
    if ((left == -Inf) & (right == Inf)) {
      p_value <- pnorm(- abs((theta - mu) / sigma))
    } else {
      if (theta != 0) {
        if (theta > 0) {
          p_value <- compute_tnorm_prob(-theta, -mu, sigma, -right, -left)
        } else {
          p_value <- compute_tnorm_prob(theta, mu, sigma, left, right)
        }
      } else {
        p_value <- NA
      }
    }
  }
  return(p_value)
}


compute_tnorm_interval <- function(theta, mu, sigma, left, right, alpha, 
                                   grid_range, grid_length) {
  grid_point <- 
    seq(from = grid_range[1] * sigma, to = grid_range[2] * sigma, length.out = grid_length)
  grid_prob <- 
    sapply(X = grid_point, FUN = function(mu_i) {
      compute_tnorm_prob(theta, mu_i, sigma, left, right) 
      },
      simplify = TRUE, 
      USE.NAMES = TRUE)
  grid_upper <- which(grid_prob >= 1 - alpha / 2)
  grid_lower <- which(grid_prob <= alpha / 2)
}


search_limit <- function(theta, mu, sigma, left, right, alpha, grid_point) {

}

compute_tnorm_prob <- function(theta, mu, sigma, left, right) {
  z_center <- (theta - mu) / sigma
  z_left <- (left - mu) / sigma
  z_right <- (right - mu) / sigma
  if (z_center <= z_left) {
    tnorm_prob <- 0
  } else if (z_center >= z_right) {
    tnorm_prob <- 1
  } else {
    tnorm_prob <- (pnorm(z_center) - pnorm(z_left)) / (pnorm(z_right) - pnorm(z_left)) 
  }
  if (is.na(tnorm_prob)) {
    if (z_left > -Inf) {
      term_left <- compute_bryc(z_left) * exp(-(z_left^2 - z_center^2) / 2)
    } else {
      term_left <- exp(z_center^2)
    }
    if (z_right < Inf) {
      term_right <- compute_bryc(z_right) * exp(-(z_right^2 - z_center^2) / 2)
    } else {
      term_right <- 0
    }
    tnorm_prob <- (term_left - compute_bryc(z_center)) / (term_left - term_right)
  }
  return(tnorm_prob)
}


compute_bryc <- function(x) {
  y <- (x^2 + 5.575192695 * x + 12.7743632) / 
    (sqrt(2 * pi) *x^3 + 14.38718147 * x^2 + 31.53531977 * x + 2 * 12.77436324)
  return(y)
}