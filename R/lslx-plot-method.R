## \code{$plot_numerical_condition()} plots the values of selected numerical conditions. ##
lslx$set("public",
         "plot_numerical_condition",
         function(condition) {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_numerical_condition' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           if (missing(condition)) {
             condition <-
               c("n_iter_out",
                 "objective_gradient_abs_max",
                 "objective_hessian_convexity")
           } else {
             if (any(!(
               condition %in% names(private$fitting$fitted_result$numerical_condition[[1]])
             ))) {
               stop("Argument `condition` contains unrecognized numerical condition.")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$numerical_condition)[condition,
                                                                                      ,
                                                                                      drop = FALSE])
           condition <-
             ifelse(
               condition == "objective_gradient_abs_max",
               "gradient",
               ifelse(
                 condition == "objective_hessian_convexity",
                 "hessian",
                 condition
               )
             )
           condition <-
             gsub(pattern = "_",
                  replacement = " ",
                  x = condition)
           
           df_for_plot$condition <- condition
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "condition",
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
           df_for_plot$delta <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
             ggplot2::geom_line() +
             ggplot2::facet_grid(condition ~ delta,
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

## \code{$plot_information_criterion()} shows how the values of information criteria vary with penalty levels. ##
lslx$set("public",
         "plot_information_criterion",
         function(criterion) {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_fit_indice' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           if (missing(criterion)) {
             criterion <- c("aic", "aic3", "caic", "bic", "abic", "hbic")
           } else {
             if (any(!(
               criterion %in% names(private$fitting$fitted_result$information_criterion[[1]])
             ))) {
               stop("Argument `criterion` contains unrecognized information criterion.")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(
               cbind,
               private$fitting$fitted_result$information_criterion
             )[criterion, , drop = FALSE])
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
           df_for_plot$delta <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
             ggplot2::geom_line(mapping = ggplot2::aes(colour = criterion)) +
             ggplot2::facet_grid(. ~ delta) +
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

## \code{$plot_fit_indice()} shows how the values of fit indices vary with penalty levels. ##
lslx$set("public",
         "plot_fit_indice",
         function(indice) {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_fit_indice' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           if (missing(indice)) {
             indice <-
               names(private$fitting$fitted_result$fit_indice[[1]])
           } else {
             if (any(!(
               indice %in% names(private$fitting$fitted_result$fit_indice[[1]])
             ))) {
               stop("Argument `indice` contains unrecognized fit indice.")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$fit_indice)[indice,
                                                                             ,
                                                                             drop = FALSE])
           df_for_plot$indice <- indice
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "indice",
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
           df_for_plot$delta <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
             ggplot2::geom_line(mapping = ggplot2::aes(colour = indice)) +
             ggplot2::facet_grid(. ~ delta) +
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

## \code{$plot_coefficient()} visualizes the solution paths of coefficients. ##
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
           df_for_plot$delta <-
             as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             ))
           ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = estimate)) +
             ggplot2::geom_line(mapping = ggplot2::aes(colour = relation)) +
             ggplot2::facet_grid(group ~ delta) +
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