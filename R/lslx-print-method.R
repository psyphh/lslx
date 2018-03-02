lslx$set("public",
         "print",
         function() {
           if (is.null(private$fitting)) {
             cat("To fit the specified model to data, please use fit-related methods.\n")
             cat("  $fit() / $fit_lasso() / $fit_mcp() / $fit_none()\n\n")
             cat("To modify the specified model, please use set-related methods.\n")
             cat("  $free_coefficient() / $penalize_coefficient() / $fix_coefficient()\n")
             cat("  $free_directed() / $penalize_directed() / $fix_directed()\n")
             cat("  $free_undirected() / $penalize_undirected() / $fix_undirected()\n")
             cat("  $free_block() / $penalize_block() / $fix_block()\n")
             cat("  $free_heterogeneity() / $penalize_heterogeneity() / $fix_heterogeneity()\n")
           } else {
             cat("To summarize the fitting results, please use $summarize().\n\n")
             cat("To plot the fitting results, please use plot-related methods.\n")
             cat("  $plot_numerical_condition() / $plot_information_criterion()\n")
             cat("  $plot_fit_indice() / $plot_coefficients()\n\n")
             cat("To obtain the test results for model or coefficients, please use test-related methods.\n")             
             cat("  $test_lr() / $test_rmsea() / $test_coefficient()\n\n")
             cat("To get a deep copy of field, please use get-related methods.\n")   
             cat("  $get_model() / $get_data() / $get_fitting()\n\n")
             cat("To extract specific quantity related to SEM, please use extract-related methods.\n")             
             cat("  $extract_specification() / $extract_saturated_mean() / $extract_saturated_cov()\n")
             cat("  $extract_saturated_moment_acov() \ $extract_penalty_level()\n")
             cat("  $extract_numerical_condition() / $extract_information_criterion()\n")
             cat("  $extract_fit_indice() / $extract_coefficient()\n")
             cat("  $extract_implied_cov() /  $extract_implied_mean()\n")
             cat("  $extract_residual_cov() / $extract_residual_mean()\n")
             cat("  $extract_coefficient_matrice() / $extract_moment_jacobian()\n")
             cat("  $extract_expected_fisher() / $extract_observed_fisher()\n")
             cat("  $extract_score_acov() / $extract_coefficient_acov()\n")
             cat("  $extract_loss_gradient() / $extract_regularizer_gradient() / $extract_objective_gradient()\n")
           }
         })
