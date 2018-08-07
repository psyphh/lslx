## \code{$plot_numerical_condition()} plots the values of selected numerical conditions. ##
lslx$set("public",
         "plot_numerical_condition",
         function(condition,
                  lambda_scale = "default",
                  mode = "default") {
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
             if (condition == "all") {
               condition <- 
                 setdiff(names(private$fitting$fitted_result$numerical_condition[[1]]),
                         c("lambda", "delta"))
             }
             if (any(!(
               condition %in% names(private$fitting$fitted_result$numerical_condition[[1]])
             ))) {
               stop("Argument `condition` contains unrecognized numerical condition.")
             }
           }
           if (lambda_scale == "default") {
             lambda_scale <- "identity"
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
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
             round(as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             )), 3)
           if (mode == "plot") {
             ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
               ggplot2::geom_line() +
               ggplot2::facet_grid(condition ~ delta,
                                   scales = "free_y",
                                   labeller = ggplot2::labeller(delta = ggplot2::label_both, 
                                                                condition = ggplot2::label_value)) +
               ggplot2::theme(
                 panel.grid.minor = ggplot2::element_line(size = .1),
                 panel.grid.major = ggplot2::element_line(size = .2)
               ) +
               ggplot2::coord_trans(x = lambda_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0("Numerical Conditions across Penalty Levels"),
                 x = "lambda",
                 y = "value"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })

## \code{$plot_information_criterion()} shows how the values of information criteria vary with penalty levels. ##
lslx$set("public",
         "plot_information_criterion",
         function(criterion,
                  lambda_scale = "default",
                  mode = "default") {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_information_criterion' method is only available for the case of 'length(lambda_grid) > 1'"
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
           if (lambda_scale == "default") {
             lambda_scale <- "identity"
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
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
             round(as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             )), 3)
           if (mode == "plot") {
             ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
               ggplot2::geom_line(mapping = ggplot2::aes(colour = criterion)) +
               ggplot2::facet_grid(. ~ delta, 
                                   labeller = ggplot2::label_both) +
               ggplot2::theme(
                 panel.grid.minor = ggplot2::element_line(size = .1),
                 panel.grid.major = ggplot2::element_line(size = .2)
               ) +
               ggplot2::coord_trans(x = lambda_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0("Values of Information Criteria across Penalty Levels"),
                 x = "lambda",
                 y = "value"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })

## \code{$plot_cv_error()} shows how the values of cv error vary with penalty levels. ##
lslx$set("public",
         "plot_cv_error",
         function(error,
                  lambda_scale = "default",
                  mode = "default") {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_cv_error()' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           if (private$fitting$control$cv_fold == 1L) {
             stop(
               "The 'plot_cv_error()' method is only available for the case of 'cv_fold > 1'."
             )
           }
           if (missing(error)) {
             error <- c("test_loss")
           } else {
             if (any(!(
               error %in% names(private$fitting$fitted_result$cv_error[[1]])
             ))) {
               stop("Argument `error` contains unrecognized cross-validation error.")
             }
           }
           if (lambda_scale == "default") {
             lambda_scale <- "identity"
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(
               cbind,
               private$fitting$fitted_result$cv_error
             )[error, , drop = FALSE])
           df_for_plot$error <- error
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "error",
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
             round(as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             )), 3)
           if (mode == "plot") {
             ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
               ggplot2::geom_line(mapping = ggplot2::aes(colour = error)) +
               ggplot2::facet_grid(. ~ delta, 
                                   labeller = ggplot2::label_both) +
               ggplot2::theme(
                 panel.grid.minor = ggplot2::element_line(size = .1),
                 panel.grid.major = ggplot2::element_line(size = .2)
               ) +
               ggplot2::coord_trans(x = lambda_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0("Values of CV Errors across Penalty Levels"),
                 x = "lambda",
                 y = "value"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })


## \code{$plot_fit_index()} shows how the values of fit indices vary with penalty levels. ##
lslx$set("public",
         "plot_fit_index",
         function(index, 
                  lambda_scale = "default",
                  mode = "default") {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_fit_index' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           if (missing(index)) {
             index <-
               names(private$fitting$fitted_result$fit_index[[1]])
           } else {
             if (any(!(
               index %in% names(private$fitting$fitted_result$fit_index[[1]])
             ))) {
               stop("Argument `index` contains unrecognized fit index.")
             }
           }
           if (lambda_scale == "default") {
             lambda_scale <- "identity"
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$fit_index)[index,
                                                                             ,
                                                                             drop = FALSE])
           df_for_plot$index <- index
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "index",
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
             round(as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             )), 3)
           if (mode == "plot") {
             ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = value)) +
               ggplot2::geom_line(mapping = ggplot2::aes(colour = index)) +
               ggplot2::facet_grid(. ~ delta, 
                                   labeller = ggplot2::label_both) +
               ggplot2::theme(
                 panel.grid.minor = ggplot2::element_line(size = .1),
                 panel.grid.major = ggplot2::element_line(size = .2)
               ) +
               ggplot2::coord_trans(x = lambda_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0("Values of Goodness-of-Fit across Penalty Levels"),
                 x = "lambda",
                 y = "value"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })

## \code{$plot_fit_indice()} shows how the values of fit indices vary with penalty levels. ##
lslx$set("public",
         "plot_fit_indice",
         function(indice,
                  lambda_scale = "default",
                  mode = "default") {
           self$plot_fit_index(indice, lambda_scale)
         })


## \code{$plot_coefficient()} visualizes the solution paths of coefficients. ##
lslx$set("public",
         "plot_coefficient",
         function(block,
                  left,
                  right,
                  both,
                  lambda_scale = "default",
                  mode = "default") {
           if (length(private$fitting$control$lambda_grid) <= 1) {
             stop(
               "The 'plot_coefficient' method is only available for the case of 'length(lambda_grid) > 1'"
             )
           }
           if (lambda_scale == "default") {
             lambda_scale <- "identity"
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
             }
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
                                   private$fitting$fitted_result$coefficient))
           df_for_plot$name <- rownames(df_for_plot)
           df_for_plot$type <- private$model$specification$type
           df_for_plot <- df_for_plot[(private$model$specification$block %in% block) &
               (private$model$specification$left %in% left) &
               (private$model$specification$right %in% right) &
               (private$model$specification$left %in% both) &
               (private$model$specification$right %in% both), , drop = FALSE]
           if (nrow(df_for_plot) == 0) {
             stop("No such type of coefficient in the specified model.")
           }
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "name",
               v.names = "estimate",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-c(ncol(df_for_plot) - 1, ncol(df_for_plot))],
               times = colnames(df_for_plot)[-c(ncol(df_for_plot) - 1, ncol(df_for_plot))],
               direction = "long"
             )
           name_split <- strsplit(x = df_for_plot$name,
                                  split = "/")
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
             round(as.numeric(sapply(
               X = penalty_level_split,
               FUN = function(x) {
                 x[4]
               }
             )), 3)
           if (mode == "plot") {
             ggplot2::ggplot(df_for_plot, ggplot2::aes(x = lambda, y = estimate)) +
               ggplot2::geom_line(mapping = ggplot2::aes(colour = relation, linetype = type),
                                  show.legend = c(colour = T, linetype = F)) +
               ggplot2::scale_linetype_manual(values = c("fixed" = "solid", 
                                                         "free" = "solid", 
                                                         "pen" = "twodash")) +
               ggplot2::facet_grid(group ~ delta, 
                                   labeller = ggplot2::label_both) +
               ggplot2::theme(
                 panel.grid.minor = ggplot2::element_line(size = .1),
                 panel.grid.major = ggplot2::element_line(size = .2)
               ) +
               ggplot2::coord_trans(x = lambda_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0(
                   "Solution Paths of Coefficients in Block ",
                   do.call(paste, as.list(block))
                 ),
                 x = "lambda",
                 y = "coefficient estimate"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })


