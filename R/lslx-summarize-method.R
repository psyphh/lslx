lslx$set("public",
         "summarize",
         function(selector,
                  standard_error = "default",
                  alpha_level = .05,
                  digit = 3,
                  exclude_improper = TRUE) {
           if (standard_error == "default") {
             if (private$fitting$control$raw) {
               standard_error <- "sandwich"
             } else {
               standard_error <- "observed_fisher"
             }
           }
           
           general_information <-
             formatC(
               x = c(
                 private$fitting$reduced_data$n_observation,
                 private$fitting$reduced_data$n_complete_observation,
                 private$fitting$reduced_data$n_missing_pattern,
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
               "number of observation",
               "number of complete observation",
               "number of missing pattern",
               "number of group",
               "number of response",
               "number of factor",
               "number of free parameter",
               "number of penalized parameter"
             )
           fitting_information <-
             formatC(
               x = c(
                 private$fitting$control$penalty_method,
                 ifelse(
                   length(private$fitting$control$lambda_grid) == 1,
                   private$fitting$control$lambda_grid,
                   paste(
                     min(private$fitting$control$lambda_grid),
                     max(private$fitting$control$lambda_grid),
                     sep = " - "
                   )
                 ),
                 ifelse(
                   length(private$fitting$control$gamma_grid) == 1,
                   private$fitting$control$gamma_grid,
                   paste(
                     min(private$fitting$control$gamma_grid),
                     max(private$fitting$control$gamma_grid),
                     sep = " - "
                   )
                 ),
                 sub(pattern = "_", replacement = " ", 
                     x = private$fitting$control$missing_method),
                 private$fitting$control$tol_out
               ),
               digits = digit,
               format = "f"
             )
           names(fitting_information) <-
             c("penalty method",
               "lambda grid",
               "gamma grid",
               "missing method",
               "tolerance for convergence")
           saturated_model_information <-
             formatC(
               x = private$fitting$supplied_result$saturated_model,
               digits = digit,
               format = "f"
             )
           names(saturated_model_information) <-
             c("loss value",
               "number of non-zero coefficient",
               "degree of freedom")
           baseline_model_information <-
             formatC(
               x = private$fitting$supplied_result$baseline_model,
               digits = digit,
               format = "f"
             )
           names(baseline_model_information) <-
             c("loss value",
               "number of non-zero coefficient",
               "degree of freedom")
           numerical_condition <-
             formatC(
               x = self$extract_numerical_condition(selector = selector,
                                                    exclude_improper = exclude_improper),
               digits = digit,
               format = "f"
             )
           names(numerical_condition) <-
             c(
               "lambda",
               "gamma",
               "objective value",
               "objective gradient abs.maximum",
               "objective hessian convexity",
               "number of iteration",
               "loss value",
               "number of non-zero coefficient",
               "degree of freedom"
             )
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
               "Akaike information criterion with penalty being 3 (aic3)",
               "consistent Akaike information criterion (caic)",
               "Bayesian information criterion (bic))",
               "adjusted Bayesian information criterion (abic)",
               "Haughton Bayesian information criterion (hbic)"
             )
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
               "comparative fit indice (cfi)",
               "non-normed fit indice (nnfi)",
               "standardized root mean of residual (srmr)"
             )
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
           
           rowname_width <-
             max(nchar(unlist(lapply(
               X = summary_list, FUN = names
             )))) + 5
           value_width <-
             max(unlist(lapply(X = summary_list, FUN = nchar))) + 3
           
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
           
           lr_test_rounded <-
             sapply(self$test_lr(selector = selector), function(i_list) {
               lr_test_single_column <-
                 formatC(i_list, digits = digit, format = "f")
               lr_test_single_column[grepl("NA", lr_test_single_column)] <-
                 "  -  "
               return(lr_test_single_column)
             })
           lr_test_rounded <- data.frame(lr_test_rounded)
           colnames(lr_test_rounded) <-
             format(c("statistic", "df", "p-value"),
                    width = 10,
                    justify = "right")
           rownames(lr_test_rounded) <-
             paste0("   ", rownames(self$test_lr(selector = selector)), "  ")
           
           cat("Likelihood Ratio Test\n")
           print(lr_test_rounded)
           cat("\n")
           
           ## printing root mean square error of approximation test
           
           rmsea_test_rounded <-
             sapply(self$test_rmsea(selector = selector), function(i_list) {
               rmsea_test_single_column <-
                 formatC(i_list, digits = digit, format = "f")
               rmsea_test_single_column[grepl("NA", rmsea_test_single_column)] <-
                 "  -  "
               return(rmsea_test_single_column)
             })
           rmsea_test_rounded <- data.frame(rmsea_test_rounded)
           colnames(rmsea_test_rounded) <-
             format(colnames(self$test_rmsea(selector = selector)),
                    width = 10,
                    justify = "right")
           rownames(rmsea_test_rounded) <-
             paste0("   ", rownames(self$test_rmsea(selector = selector)), "  ")
           
           cat("Root Mean Square Error of Approximation Test\n")
           print(rmsea_test_rounded)
           cat("\n")
           
           
           
           ## printing coefficient
           
           coefficient <- self$test_coefficient(
             selector = selector,
             standard_error = standard_error,
             exclude_improper = exclude_improper
           )
           
           relation_as_groupname <-
             format(
               private$model$specification$relation,
               width = max(nchar(private$model$specification$relation)),
               justify = "right"
             )
           block_levels <-
             c("Factor Loading",
               "Regression",
               "Covariance",
               "Variance",
               "Intercept")
           
           coefficient_rounded <-
             lapply(coefficient, function(i_list) {
               data_single_column <- formatC(i_list, digits = digit, format = "f")
               data_single_column[grepl("NA", data_single_column)] <-
                 "  -  "
               return(data_single_column)
             })
           
           coefficient_rounded$block <-
             private$model$specification$block
           coefficient_rounded$block_type <-
             rep(NA, nrow(private$model$specification))
           coefficient_rounded$type <-
             format(private$model$specification$type,
                    width = 6,
                    justify = "right")
           coefficient_rounded$group <-
             private$model$specification$group
           coefficient_rounded$left <-
             private$model$specification$left
           coefficient_rounded$right <-
             private$model$specification$right
           
           coefficient_rounded$block_type[coefficient_rounded$block == "y<-1"] <-
             "Intercept"
           coefficient_rounded$block_type[coefficient_rounded$block == "f<-1"] <-
             "Intercept"
           coefficient_rounded$block_type[coefficient_rounded$block == "y<-f"] <-
             "Factor Loading"
           coefficient_rounded$block_type[coefficient_rounded$block == "y<-y"] <-
             "Regression"
           coefficient_rounded$block_type[coefficient_rounded$block == "f<-y"] <-
             "Regression"
           coefficient_rounded$block_type[coefficient_rounded$block == "f<-f"] <-
             "Regression"
           coefficient_rounded$block_type[coefficient_rounded$block == "y<->f"] <-
             "Covariance"
           coefficient_rounded$block_type[coefficient_rounded$block == "f<->y"] <-
             "Covariance"
           coefficient_rounded$block_type[coefficient_rounded$block == "y<->y"] <-
             "Covariance"
           coefficient_rounded$block_type[coefficient_rounded$block == "f<->f"] <-
             "Covariance"
           coefficient_rounded$block_type[coefficient_rounded$left == coefficient_rounded$right] <-
             "Variance"
           
           coefficient_rounded <-
             data.frame(coefficient_rounded[c(1:10)])
           
           ## print by different groups
           if (!is.na(private$model$reference_group)) {
             reference_group_order <-
               which(private$model$name_group %in% private$model$reference_group)
             group_by_order <- c(reference_group_order,
                                 c(1:(length(
                                   private$model$name_group
                                 )))[!(c(1:(length(
                                   private$model$name_group
                                 ))) %in% reference_group_order)])
           } else {
             group_by_order <- 1:length(private$model$name_group)
           }
           
           for (i_group in group_by_order) {
             selected_rows <-
               coefficient_rounded$group == private$model$name_group[i_group]
             
             data_single_group <-
               coefficient_rounded[selected_rows, ]
             
             rownames(data_single_group) <-
               paste0(relation_as_groupname[selected_rows], "  ")
             
             colnames(data_single_group) <- c(
               "estimate",
               "std.error",
               "z-value",
               "p-value",
               "lower",
               "upper",
               "block",
               "block_type",
               "type",
               "group"
             )
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
                 cat(" ", i_block_type, "\n")
                 dta_single_group <-
                   data_single_group[data_single_group$block_type == i_block_type, c(9, 1:6)]
                 colnames(dta_single_group) <-
                   format(colnames(dta_single_group),
                          width = 8,
                          justify = "right")
                 print(dta_single_group)
                 cat("\n")
               }
             }
           }
         })
