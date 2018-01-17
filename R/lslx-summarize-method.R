## \code{$summarize()} prints a summary for the fitting result under the given selector. ##
lslx$set("public",
         "summarize",
         function(selector,
                  standard_error = "default",
                  alpha_level = .05,
                  digit = 3,
                  interval = TRUE,
                  simplify = FALSE,
                  exclude_improper = TRUE) {
           setting <- list(
             general_information = TRUE,
             fitting_information = TRUE,
             saturated_model_information = TRUE,
             baseline_model_information = TRUE,
             numerical_condition = TRUE,
             information_criterion = TRUE,
             fit_indices = TRUE,
             likelihood_ratio_test = TRUE,
             root_mean_square_error_of_approximation_test = TRUE,
             coefficient_test = TRUE
           )
           
           coefficient_test_setting <- list(
             type = TRUE,
             estimate = TRUE,
             standard_error = TRUE,
             z_value = TRUE,
             p_value = TRUE,
             lower = TRUE,
             upper = TRUE,
             block = FALSE,
             #block and group default don't print (as columns)
             block_type = TRUE,
             # the block_type is an index which should be removed at last
             group = FALSE
           ) #block and group default don't print (as columns)
           
           if (simplify) {
             setting$fitting_information <- FALSE
             setting$saturated_model_information <- FALSE
             setting$baseline_model_information <- FALSE
             setting$information_criterion <- FALSE
             setting$fit_indices <- FALSE
             setting$likelihood_ratio_test <- FALSE
             setting$root_mean_square_error_of_approximation_test <-
               FALSE
             setting$coefficient_test <- FALSE
           }
           
           if (!interval) {
             coefficient_test_setting$upper <- FALSE
             coefficient_test_setting$lower <- FALSE
           }

           if (standard_error == "default") {
             if (private$fitting$control$response) {
               standard_error <- "sandwich"
             } else {
               standard_error <- "observed_fisher"
             }
           }
           
           ##generating output informations
           if (setting$general_information) {
             general_information <-
               formatC(
                 x = c(
                   private$fitting$reduced_data$n_observation,
                   private$fitting$reduced_data$n_complete_observation,
                   ifelse(
                     private$fitting$reduced_data$n_missing_pattern == 1,
                     "none",
                     private$fitting$reduced_data$n_missing_pattern
                   ),
                   private$fitting$reduced_model$n_group,
                   private$fitting$reduced_model$n_response,
                   private$fitting$reduced_model$n_factor,
                   sum(private$fitting$reduced_model$theta_is_free),
                   sum(private$fitting$reduced_model$theta_is_pen)
                 ),
                 digits = digit,
                 format = "f"
               )
             names(general_information) <-
               c(
                 "number of observations",
                 "number of complete observations",
                 "number of missing patterns",
                 "number of groups",
                 "number of responses",
                 "number of factors",
                 "number of free coefficients",
                 "number of penalized coefficients"
               )
           } else {
             general_information <- NULL
           }
           
           if (setting$fitting_information) {
             fitting_information <-
               formatC(
                 x = c(
                   private$fitting$control$penalty_method,
                   ifelse(
                     private$fitting$control$penalty_method == "none",
                     "none",
                     ifelse(
                       length(private$fitting$control$lambda_grid) == 1,
                       private$fitting$control$lambda_grid,
                       paste(
                         min(private$fitting$control$lambda_grid),
                         max(private$fitting$control$lambda_grid),
                         sep = " - "
                       )
                     )
                   ),
                   ifelse(
                     private$fitting$control$penalty_method != "mcp",
                     "none",
                     ifelse(
                       length(private$fitting$control$delta_grid) == 1,
                       private$fitting$control$delta_grid,
                       paste(
                         min(private$fitting$control$delta_grid),
                         max(private$fitting$control$delta_grid),
                         sep = " - "
                       )
                     )
                   ),
                   private$fitting$control$algorithm,
                   ifelse(
                     private$fitting$reduced_data$n_missing_pattern == 1,
                     "none",
                     sub(
                       pattern = "_",
                       replacement = " ",
                       x = private$fitting$control$missing_method
                     )
                   ),
                   private$fitting$control$tol_out
                 ),
                 digits = digit,
                 format = "f"
               )
             names(fitting_information) <-
               c(
                 "penalty method",
                 "lambda grid",
                 "delta grid",
                 "algorithm",
                 "missing method",
                 "tolerance for convergence"
               )
           } else {
             fitting_information <- NULL
           }
           
           if (setting$saturated_model_information) {
             saturated_model_information <-
               formatC(
                 x = private$fitting$supplied_result$saturated_model,
                 digits = digit,
                 format = "f"
               )
             names(saturated_model_information) <-
               c("loss value",
                 "number of non-zero coefficients",
                 "degree of freedom")
           } else {
             saturated_model_information <- NULL
           }
           
           if (setting$baseline_model_information) {
             baseline_model_information <-
               formatC(
                 x = private$fitting$supplied_result$baseline_model,
                 digits = digit,
                 format = "f"
               )
             names(baseline_model_information) <-
               c("loss value",
                 "number of non-zero coefficients",
                 "degree of freedom")
           } else {
             baseline_model_information <- NULL
           }
           
           if (setting$numerical_condition) {
             numerical_condition <-
               formatC(
                 x = self$extract_numerical_condition(selector = selector,
                                                      exclude_improper = exclude_improper),
                 digits = digit,
                 format = "f"
               )
             numerical_condition[["lambda"]] <-
               ifelse(private$fitting$control$penalty_method == "none",
                      "none",
                      numerical_condition[["lambda"]])
             numerical_condition[["delta"]] <-
               ifelse(private$fitting$control$penalty_method != "mcp",
                      "none",
                      numerical_condition[["delta"]])
             names(numerical_condition) <-
               c(
                 "lambda",
                 "delta",
                 "objective value",
                 "objective gradient absolute maximum",
                 "objective Hessian convexity",
                 "number of iterations",
                 "loss value",
                 "number of non-zero coefficients",
                 "degree of freedom",
                 "robust degree of freedom",
                 "scaling factor"
               )
           } else {
             numerical_condition <- NULL
           }
           
           if (setting$information_criterion) {
             information_criterion <-
               formatC(
                 x = self$extract_information_criterion(selector = selector,
                                                        exclude_improper = exclude_improper),
                 digits = digit,
                 format = "f"
               )
             names(information_criterion) <-
               c(
                 "Akaike information criterion (aic)",
                 "Akaike information criterion with penalty 3 (aic3)",
                 "consistent Akaike information criterion (caic)",
                 "Bayesian information criterion (bic)",
                 "adjusted Bayesian information criterion (abic)",
                 "Haughton Bayesian information criterion (hbic)",
                 "robust Akaike information criterion (raic)",
                 "robust Akaike information criterion with penalty 3 (raic3)",
                 "robust consistent Akaike information criterion (rcaic)",
                 "robust Bayesian information criterion (rbic)",
                 "robust adjusted Bayesian information criterion (rabic)",
                 "robust Haughton Bayesian information criterion (rhbic)"
               )
           } else {
             information_criterion <- NULL
           }
           
           if (setting$fit_indice) {
             fit_indice <-
               formatC(
                 x = self$extract_fit_indice(selector = selector,
                                             exclude_improper = exclude_improper),
                 digits = digit,
                 format = "f"
               )
             names(fit_indice) <-
               c(
                 "root mean square error of approximation (rmsea)",
                 "comparative fit index (cfi)",
                 "non-normed fit index (nnfi)",
                 "standardized root mean of residual (srmr)"
               )
           } else {
             fit_indice <- NULL
           }
           
           if (setting$likelihood_ratio_test) {
             lr_test <-
               self$test_lr(selector = selector,
                            exclude_improper = exclude_improper)
             lr_test_rounded <-
               data.frame(sapply(
                 X = lr_test,
                 FUN = function(lr_test_i) {
                   lr_test_rounded_i <-
                     formatC(lr_test_i, digits = digit, format = "f")
                   lr_test_rounded_i[grepl("NA", lr_test_rounded_i)] <-
                     "  -  "
                   return(lr_test_rounded_i)
                 }
               ))
             colnames(lr_test_rounded) <-
               format(c("statistic", "df", "p-value"),
                      width = 10,
                      justify = "right")
             rownames(lr_test_rounded) <-
               paste0("   ", rownames(lr_test), "  ")
           }
           
           if (setting$root_mean_square_error_of_approximation_test) {
             rmsea_test <-
               self$test_rmsea(
                 selector = selector,
                 alpha_level = alpha_level,
                 exclude_improper = exclude_improper
               )
             rmsea_test_rounded <-
               data.frame(sapply(
                 X = rmsea_test,
                 FUN = function(rmsea_test_i) {
                   rmsea_test_rounded_i <-
                     formatC(rmsea_test_i, digits = digit, format = "f")
                   rmsea_test_rounded_i[grepl("NA", rmsea_test_rounded_i)] <-
                     "  -  "
                   return(rmsea_test_rounded_i)
                 }
               ))
             colnames(rmsea_test_rounded) <-
               format(colnames(rmsea_test),
                      width = 10,
                      justify = "right")
             rownames(rmsea_test_rounded) <-
               paste0("   ", rownames(rmsea_test), "  ")
           }

           summary_list <-
             list(
               general_information,
               fitting_information,
               saturated_model_information,
               baseline_model_information,
               numerical_condition,
               information_criterion,
               fit_indice
             )
           names(summary_list) <- c(
             "General Information",
             "Fitting Information",
             "Saturated Model Information",
             "Baseline Model Information",
             "Numerical Condition",
             "Information Criteria",
             "Fit Indices"
           )
           
           summary_list <-
             summary_list[!sapply(summary_list, is.null)]
           
           rowname_width <-
             max(nchar(unlist(lapply(
               X = summary_list, FUN = names
             )))) + 5
           value_width <-
             max(unlist(lapply(X = summary_list, FUN = nchar))) + 1
           
           ## printing summary list
           for (name_i in names(summary_list)) {
             cat(name_i)
             summary_list_i <-
               as.data.frame(summary_list[[name_i]])
             colnames(summary_list_i) <- NULL
             rownames(summary_list_i) <-
               format(paste("  ", rownames(summary_list_i)),
                      width = rowname_width,
                      justify = "left")
             print(format(summary_list_i, width = value_width, justify = "right"))
             cat("\n")
           }
           
           
           
           ## printing likelihood ratio test
           if (setting$likelihood_ratio_test) {
             cat("Likelihood Ratio Test\n")
             print(lr_test_rounded)
             cat("\n")
           }
           
           ## printing root mean square error of approximation test
           if (setting$root_mean_square_error_of_approximation_test) {
             cat("Root Mean Square Error of Approximation Test\n")
             print(rmsea_test_rounded)
             cat("\n")
           }
           
           ## generating coefficients
           if (setting$coefficient_test) {
             coefficient_test <-
               self$test_coefficient(
                 selector = selector,
                 standard_error = standard_error,
                 alpha_level = alpha_level,
                 exclude_improper = exclude_improper
               )
             relation_as_groupname <-
               format(
                 private$model$specification$relation,
                 width = max(nchar(
                   private$model$specification$relation
                 )),
                 justify = "right"
               )
             block_levels <-
               c("Factor Loading",
                 "Regression",
                 "Covariance",
                 "Variance",
                 "Intercept")
             coefficient_test_rounded <-
               lapply(
                 X = coefficient_test,
                 FUN = function(coefficient_test_i) {
                   coefficient_test_rounded_i <-
                     formatC(coefficient_test_i,
                             digits = digit,
                             format = "f")
                   coefficient_test_rounded_i[grepl("NA", coefficient_test_rounded_i)] <-
                     "  -  "
                   return(coefficient_test_rounded_i)
                 }
               )
             
             coefficient_test_rounded <-
               c(coefficient_test_rounded, private$model$specification[c("block",
                                                                         "group",
                                                                         "left",
                                                                         "right")])
             
             coefficient_test_rounded$block_type <-
               rep(NA, nrow(private$model$specification))
             coefficient_test_rounded$type <-
               format(private$model$specification$type,
                      width = 6,
                      justify = "right")
             
             
             coefficient_test_rounded$block_type[grepl("(y|f)<-1", coefficient_test_rounded$block)] <-
               "Intercept"
             coefficient_test_rounded$block_type[grepl("(y|f)<-(y|f)", coefficient_test_rounded$block)] <-
               "Regression"
             coefficient_test_rounded$block_type[coefficient_test_rounded$block == "y<-f"] <-
               "Factor Loading"
             coefficient_test_rounded$block_type[grepl("(y|f)<->(y|f)", coefficient_test_rounded$block)] <-
               "Covariance"
             coefficient_test_rounded$block_type[coefficient_test_rounded$left == coefficient_test_rounded$right] <-
               "Variance"
             
             coefficient_test_rounded <-
               data.frame(coefficient_test_rounded[c(
                 "type",
                 "estimate",
                 "standard_error",
                 "z_value",
                 "p_value",
                 "lower",
                 "upper",
                 "block",
                 "block_type",
                 "group"
               )])
             
             ## print by different groups
             if (!is.na(private$model$reference_group)) {
               reference_group_order <-
                 which(private$model$name_group %in% private$model$reference_group)
               group_by_order <-
                 c(reference_group_order,
                   c(1:(length(
                     private$model$name_group
                   )))[!(c(1:(length(
                     private$model$name_group
                   ))) %in% reference_group_order)])
             } else {
               group_by_order <- 1:length(private$model$name_group)
             }
             
             for (i_group in group_by_order) {
               idc_group <-
                 coefficient_test_rounded$group == private$model$name_group[i_group]
               data_single_group <-
                 coefficient_test_rounded[idc_group,]
               rownames(data_single_group) <-
                 paste0(relation_as_groupname[idc_group], "  ")
               colnames(data_single_group) <-
                 c(
                   "type",
                   "estimate",
                   "std.error",
                   "z-value",
                   "p-value",
                   "lower",
                   "upper",
                   "block",
                   "block_type",
                   "group"
                 )
               
               data_single_group <-
                 data_single_group[unlist(coefficient_test_setting)]
               
               
               if (length(private$model$name_group) != 1) {
                 cat(
                   paste0(
                     "Coefficient Test (Group = \"",
                     private$model$name_group[[i_group]],
                     "\"",
                     ", Standard Error = \"",
                     standard_error,
                     "\"",
                     ", Alpha Level = ",
                     alpha_level,
                     ")\n"
                   )
                 )
               } else {
                 cat(
                   paste0(
                     "Coefficient Test",
                     " (Standard Error = \"",
                     standard_error,
                     "\"",
                     ", Alpha Level = ",
                     alpha_level,
                     ")\n"
                   )
                 )
               }
               
               ## print by block types
               for (i_block_type in block_levels) {
                 if (sum(data_single_group$block_type == i_block_type) > 0L) {
                   cat(" ", i_block_type)
                   # if 'single group' or 'reference group not specified', print nothing.
                   if ((length(group_by_order) == 1) |
                       is.na(private$model$reference_group)) {
                     cat("\n")
                   } else if (i_group == group_by_order[1]) {
                     cat(" (reference component)\n")
                   } else {
                     cat(" (increment component)\n")
                   }
                   data_single_group_block <-
                     data_single_group[(data_single_group$block_type == i_block_type),!(colnames(data_single_group) %in% c("block_type"))]
                   colnames(data_single_group_block) <-
                     paste0(" ", colnames(data_single_group_block))
                   print(data_single_group_block)
                   cat("\n")
                 }
               }
             }
           }
         })
