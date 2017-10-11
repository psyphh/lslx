lslx$set("public",
         "plot_numerical_condition",
         function() {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_numerical_condition' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }

           criterion <-
             c("n_iter_out",
               "objective_gradient_abs_max",
               "objective_hessian_convexity")
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$numerical_condition)[criterion,
                                                                                      ,
                                                                                      drop = FALSE])
           criterion <-
             gsub(
               pattern = "_",
               replacement = " ",
               x = rownames(df_for_plot)
             )
           df_for_plot$criterion <- criterion
           
           
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "criterion",
               v.names = "value",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-ncol(df_for_plot)],
               times = colnames(df_for_plot)[-ncol(df_for_plot)],
               direction = "long"
             )
           
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|/")
           
           df_for_plot$penalty_level <- NULL
           
           df_for_plot$lambda <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[2]
               }
             ))
           
           df_for_plot$gamma <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
             ggplot2::geom_line() +
             ggplot2::facet_grid(criterion ~ gamma,
                                 scales = "free_y") +
             ggplot2::theme(
               panel.grid.minor = ggplot2::element_line(size = .1),
               panel.grid.major = ggplot2::element_line(size = .2)
             )  +
             ggplot2::labs(
               title = paste0("Numerical Conditions across Penalty Levels"),
               x = "lambda",
               y = "value"
             ) +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
         })

lslx$set("public",
         "plot_information_criterion",
         function() {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_fit_indice' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           selector <-
             names(private$fitting$fitted_result$information_criterion[[1]])
           
           df_for_plot <-
             as.data.frame(
               do.call(
                 cbind,
                 private$fitting$fitted_result$information_criterion)[
                   selector, , drop = FALSE])
           
           df_for_plot$selector <- selector
           
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "selector",
               v.names = "value",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-ncol(df_for_plot)],
               times = colnames(df_for_plot)[-ncol(df_for_plot)],
               direction = "long"
             )
           
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|/")
           
           df_for_plot$penalty_level <- NULL
           
           df_for_plot$lambda <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[2]
               }
             ))
           
           df_for_plot$gamma <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
             ggplot2::geom_line(mapping = ggplot2::aes(colour = selector)) +
             ggplot2::facet_grid(. ~ gamma, labeller = ggplot2::label_both) +
             ggplot2::theme(
               panel.grid.minor = ggplot2::element_line(size = .1),
               panel.grid.major = ggplot2::element_line(size = .2)
             )  +
             ggplot2::labs(
               title = paste0("Values of Information Criteria across Penalty Levels"),
               x = "lambda",
               y = "value"
             ) +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
         })



lslx$set("public",
         "plot_fit_indice",
         function() {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_fit_indice' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           selector <-
             names(private$fitting$fitted_result$fit_indice[[1]])
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$fit_indice)[selector,
                                                                                  ,
                                                                                  drop = FALSE])
           df_for_plot$selector <- selector
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "selector",
               v.names = "value",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-ncol(df_for_plot)],
               times = colnames(df_for_plot)[-ncol(df_for_plot)],
               direction = "long"
             )
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|/")
           df_for_plot$penalty_level <- NULL
           df_for_plot$lambda <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[2]
               }
             ))
           df_for_plot$gamma <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
             ggplot2::geom_line(mapping = ggplot2::aes(colour = selector)) +
             ggplot2::facet_grid(. ~ gamma, labeller = ggplot2::label_both) +
             ggplot2::theme(
               panel.grid.minor = ggplot2::element_line(size = .1),
               panel.grid.major = ggplot2::element_line(size = .2)
             )  +
             ggplot2::labs(
               title = paste0("Values of Goodness-of-Fit across Penalty Levels"),
               x = "lambda",
               y = "value"
             ) +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
         })


lslx$set("public",
         "plot_coefficient",
         function(block,
                  left,
                  right,
                  both) {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_coefficient' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           if (missing(block)) {
             block <- unique(private$model$specification$block)
           }
           if (missing(left)) {
             left <- c(private$model$name_eta, "1")
           }
           if (missing(right)) {
             right <- c(private$model$name_eta, "1")
           }
           if (missing(both)) {
             both <- c(private$model$name_eta, "1")
           }
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$coefficient)[(private$model$specification$block %in% block) &
                                                                                (private$model$specification$left %in% left) &
                                                                                (private$model$specification$right %in% right) &
                                                                                (private$model$specification$left %in% both) &
                                                                                (private$model$specification$right %in% both),
                                                                              ,
                                                                              drop = FALSE])
           
           if (nrow(df_for_plot) == 0) {
             sleftp("No such type of coefficient in the specified model.")
           }
           df_for_plot$name <- rownames(df_for_plot)
           
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "name",
               v.names = "estimate",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-ncol(df_for_plot)],
               times = colnames(df_for_plot)[-ncol(df_for_plot)],
               direction = "long"
             )
           
           name_split <- strsplit(x = df_for_plot$name,
                                  split = "\\|")
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|/")
           
           df_for_plot$name <- NULL
           df_for_plot$penalty_level <- NULL
           
           df_for_plot$relation <-
             sapply(
               X = name_split,
               FUN = function(x) {
                 x[1]
               }
             )
           
           df_for_plot$group <-
             sapply(
               X = name_split,
               FUN = function(x) {
                 x[2]
               }
             )
           
           df_for_plot$lambda <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[2]
               }
             ))
           
           df_for_plot$gamma <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = estimate)) +
             ggplot2::geom_line(mapping = ggplot2::aes(colour = relation)) +
             ggplot2::facet_grid(group ~ gamma, labeller = ggplot2::label_both) +
             ggplot2::theme(
               panel.grid.minor = ggplot2::element_line(size = .1),
               panel.grid.major = ggplot2::element_line(size = .2)
             )  +
             ggplot2::labs(
               title = paste0(
                 "Solution Paths of Coefficients in Block ",
                 do.call(paste, as.list(block))
               ),
               x = "lambda",
               y = "coefficient estimate"
             ) +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
         })