lslx$set("public",
         "summarize",
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
           
           
           fitted_coefficient <-
             self$extract_coefficient(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )[[1]]
           
           
           general_information <- list()
           general_information$sample_size <- private$fitting$reduced_data$total_sample_size 
           general_information$n_group <- private$fitting$reduced_model$n_group
           general_information$n_response <- private$fitting$reduced_model$n_response
           general_information$n_factor <- private$fitting$reduced_model$n_factor
           general_information$n_free_parameter <- sum(private$fitting$reduced_model$theta_is_free)
           general_information$n_penalized_parameter <- sum(private$fitting$reduced_model$theta_is_pen)
           
           fitting_information <- list()
           fitting_information$penalty_method <- private$fitting$control$penalty_method
           fitting_information$lambda_grid <- private$fitting$control$lambda_grid
           fitting_information$gamma_grid <- private$fitting$control$gamma_grid
           
           saturated_model_information <- list()           
           saturated_model_information$loss_value <- 0
           saturated_model_information$n_nonzero_coefficient <- 
             private$fitting$reduced_model$n_moment * private$fitting$reduced_model$n_group
           
           baseline_model_information <- list()           
           baseline_model_information$loss_value <- private$fitting$reduced_model$loss_value_baseline
           baseline_model_information$n_nonzero_coefficient <- 
             2 * private$fitting$reduced_model$n_response * private$fitting$reduced_model$n_group
           
           selected_model_information <- list()           
           selected_model_information$selector <- selector
           selected_model_information$lambda <- fitted_numerical_condition[["lambda"]]
           selected_model_information$gamma <- fitted_numerical_condition[["gamma"]]
           selected_model_information$n_iter_out <- fitted_numerical_condition[["n_iter_out"]]
           selected_model_information$objective_value <- fitted_numerical_condition[["objective_value"]]
           selected_model_information$loss_value <- fitted_goodness_of_fit[["loss"]]
           selected_model_information$n_nonzero_coefficient <- fitted_numerical_condition[["n_nonzero_coefficient"]]
           
           information_criteria <- list()
           information_criteria$aic <- fitted_goodness_of_fit[["aic"]]
           information_criteria$aic3 <- fitted_goodness_of_fit[["aic3"]]
           information_criteria$caic <- fitted_goodness_of_fit[["caic"]]
           information_criteria$bic <- fitted_goodness_of_fit[["bic"]]
           information_criteria$abic <- fitted_goodness_of_fit[["abic"]]
           information_criteria$hbic <- fitted_goodness_of_fit[["hbic"]]
           
           likelihood_ratio <-
             self$test_likelihood_ratio(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           
           root_mean_square_error_of_approximation <-
             self$test_rmsea(
               selector = selector,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           
           comparative_fit_indice <- list()
           comparative_fit_indice$cfi <- fitted_goodness_of_fit[["cfi"]]
           
           non_normed_fit_indice <- list()
           non_normed_fit_indice$nnfi <- fitted_goodness_of_fit[["nnfi"]]
           
           standardized_root_mean_of_residual <- list() 
           standardized_root_mean_of_residual$srmr <- fitted_goodness_of_fit[["srmr"]]
           
           coefficient <-
             self$test_coefficient(
               selector = selector,
               standard_error = standard_error,
               exclude_nonconvergence = exclude_nonconvergence,
               exclude_nonconvexity = exclude_nonconvexity,
               verbose = verbose
             )
           list_summary <- list(general_information = general_information,
                                fitting_information = fitting_information,
                                saturated_model_information = saturated_model_information,
                                baseline_model_information = baseline_model_information,
                                selected_model_information = selected_model_information,
                                information_criteria = information_criteria,
                                likelihood_ratio = likelihood_ratio,
                                root_mean_square_error_of_approximation = root_mean_square_error_of_approximation,
                                comparative_fit_indice = comparative_fit_indice,
                                non_normed_fit_indice = non_normed_fit_indice,
                                standardized_root_mean_of_residual = standardized_root_mean_of_residual,
                                coefficient = coefficient)
           
           
           digit = 3
           
           
           relation_as_groupname <- format(private$model$specification$relation,
                                           width = max(nchar(private$model$specification$relation)),
                                           justify = "right")
           block_levels <- c("Factor Loading","Regression","Covariance","Variance","Intercept")
           
           fitting_information$lambda_grid <- paste(min(fitting_information$lambda_grid),
                                                    max(fitting_information$lambda_grid),
                                                    sep = " - ")
           fitting_information$gamma_grid <- paste(min(fitting_information$gamma_grid),
                                                   max(fitting_information$gamma_grid),
                                                   sep = " - ")
           
           
           information_names <- c(gen_inf = "General Information",
                                  fit_inf = "Fitting Information",
                                  sat_mod_inf = "Saturated Model Information",
                                  bas_mod_inf = "Baseline Model Information",
                                  sel_mod_inf = "Selected Model Information",
                                  info_cri = "Information Crirteria",
                                  lik_rat = "Likelihood Ratio",
                                  rmsea = "Root Mean Square Error of Approximation",
                                  cfi = "Comparative Fit Indice",
                                  nnfi = "Non Normed Fit Indice",
                                  srmr = "Standardized Root Mean of Residual")
           
           information_data <- list(gen_inf = unlist(general_information),
                                    fit_inf = unlist(fitting_information),
                                    sat_mod_inf = formatC(unlist(saturated_model_information),
                                                          digits = digit,
                                                          format = "f"),
                                    bas_mod_inf = formatC(unlist(baseline_model_information),
                                                          digits = digit,
                                                          format = "f"),
                                    sel_mod_inf = c(selector = selected_model_information[[1]],
                                                    formatC(unlist(selected_model_information[2:7]),
                                                            digits = digit,
                                                            format = "f")),
                                    info_cri = formatC(unlist(information_criteria),
                                                       digits = digit,
                                                       format = "f"),
                                    lik_rat = formatC(unlist(likelihood_ratio),
                                                      digits = digit,
                                                      format = "f"),
                                    rmsea = formatC(unlist(root_mean_square_error_of_approximation),
                                                    digits = digit,
                                                    format = "f"),
                                    cfi = formatC(unlist(comparative_fit_indice),
                                                  digits = digit,
                                                  format = "f"),
                                    nnfi = formatC(unlist(non_normed_fit_indice),
                                                   digits = digit,
                                                   format = "f"),
                                    srmr = formatC(unlist(standardized_root_mean_of_residual),
                                                   digits = digit,
                                                   format = "f"))
           
           information_row_names <- list(gen_inf = c("sample size",
                                                     "number of groups",
                                                     "number of response",
                                                     "number of factor",
                                                     "number of free parameter",
                                                     "number of penalized parameter"),
                                         fit_inf = c("penalty method",
                                                     "lambda grid",
                                                     "gamma grid"),
                                         sat_mod_inf = c("loss value",
                                                         "number of non-zero coefficient"),
                                         bas_mod_inf = c("loss value",
                                                         "number of non-zero coefficient"),
                                         sel_mod_inf = c("selector",
                                                         "lambda",
                                                         "gamma",
                                                         "number of iteration",
                                                         "objective value",
                                                         "loss value",
                                                         "number of non-zero coefficient"),
                                         info_cri =  c("aic (Akaike, 1974)",
                                                       "aic3 (Sclove, 1987)",
                                                       "caic (Bozdogan, 1987)",
                                                       "bic (Schwarz, 1978)",
                                                       "abic (Yang, 2006)",
                                                       "hbic (Haughton, 1997)"),
                                         lik_rat = c("statistic",
                                                     "degree of freedom",
                                                     "p-value"),
                                         rmsea = c("rmsea",
                                                   "lower limit",
                                                   "upper limit"),
                                         cfi = c("cfi"),
                                         nnfi = c("nnfi"),
                                         srmr = c("srmr")
           )
           
           rowname_width <- max(nchar(information_names))+5
           value_width <- max(unlist(lapply(information_data,nchar)))+3
           
           for (i_information in names(information_names)){
             cat(information_names[i_information])
             values <- as.data.frame(information_data[i_information][[1]])
             colnames(values) <- NULL
             rownames(values) <- format(paste("  ",information_row_names[i_information][[1]]),
                                        width = rowname_width,
                                        justify = "left")               
             print(format(values, width = value_width, justify = "right"))
             cat("\n")
             
           }
           
           rounded_coe <- lapply(coefficient,function(i_list) {
             single_col_dta <- formatC(i_list,digits = 3,format = "f")
             single_col_dta[grepl("NA",single_col_dta)] <- "-"
             return(single_col_dta)
           }
           )
           
           rounded_coe$block <- private$model$specification$block
           rounded_coe$block_type <- rep(NA,nrow(private$model$specification))
           rounded_coe$type <- format(private$model$specification$type,width = 6, justify = "right")
           rounded_coe$left <- private$model$specification$left
           rounded_coe$right <- private$model$specification$right
           
           
           
           rounded_coe$block_type[rounded_coe$block=="y<-1"] <- "Intercept"
           rounded_coe$block_type[rounded_coe$block=="f<-1"] <- "Intercept"
           rounded_coe$block_type[rounded_coe$block=="y<-f"] <- "Factor Loading"
           rounded_coe$block_type[rounded_coe$block=="y<-y"] <- "Regression"
           rounded_coe$block_type[rounded_coe$block=="f<-y"] <- "Regression"
           rounded_coe$block_type[rounded_coe$block=="f<-f"] <- "Regression"
           rounded_coe$block_type[rounded_coe$block=="y<->f"] <- "Covariance"
           rounded_coe$block_type[rounded_coe$block=="f<->y"] <- "Covariance"
           rounded_coe$block_type[rounded_coe$block=="y<->y"] <- "Covariance"
           rounded_coe$block_type[rounded_coe$block=="f<->f"] <- "Covariance"
           rounded_coe$block_type[rounded_coe$left==rounded_coe$right] <- "Variance"
           
           rounded_coe <- data.frame(rounded_coe[c(1:7)])
           
           
           ## print by different groups
           
           if (!is.na(private$model$reference_group)){
             group_by_order <- c(grep(private$model$reference_group,private$model$name_group),
                                 c(1:(length(private$model$name_group)))[!(c(1:(length(private$model$name_group)) %in%
                                                                               c(grep(private$model$reference_group,private$model$name_group))))])
           } else {
             group_by_order <- 1:length(private$model$name_group)
           }
           
           for (i_group in group_by_order){
             
             single_group_dta <- rounded_coe[grepl(private$model$name_group[[i_group]],
                                                   rownames(rounded_coe)),]
             
             rownames(single_group_dta) <- relation_as_groupname[grepl(private$model$name_group[[i_group]],
                                                                       rownames(rounded_coe))]
             
             colnames(single_group_dta) <- c("estimate",
                                             "std.error",
                                             "z-value",
                                             "p-value",
                                             "block",
                                             "block_type",
                                             "type")
             
             if (length(private$model$name_group)!=1) {
               cat(paste0(
                 "Coefficient Estimate (Group = \"",
                 private$model$name_group[[i_group]],
                 "\")\n"))
             } else {
               cat("Coefficient Estimate",
                   "\n")
             }
             
             
             
             ## print by block types
             
             for (i_block_type in block_levels){
               if(
                 sum(single_group_dta$block_type == i_block_type)>0L
               ) {
                 cat(" ",i_block_type,"\n")
                 print(single_group_dta[single_group_dta$block_type == i_block_type,c(7,1:4)])
                 cat("\n")
               }
             }
           }
           
           
           
           #return(list_summary)
         })
