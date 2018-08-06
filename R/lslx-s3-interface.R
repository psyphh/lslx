#' S3 interface for semi-confirmatory SEM via PL
#' 
#' \code{plsem()} is an \code{S3} interface for obaining a fitted \code{lslx} object. 
#' 
#' @return A fitted \code{lslx} object
#' @param model A \code{character} with length one to represent the model specification. 
#' @param data A \code{data.frame} of raw data.
#' @param penalty_method A \code{character} to specify the penalty method.
#' @param lambda_grid A non-negative \code{numeric} for specifying penalty levels for both \code{"lasso"} and \code{"mcp"}.
#' @param delta_grid A non-negative \code{numeric} for specifying convexity  levels for \code{"mcp"}.
#' @param ... Other arguments. For details, please see the documentation of \code{lslx}.
#' @examples 
#' ## Semi-Confirmatory Factor Analysis with lavaan Style ##
#' # specify a factor analysis model with lavaan style
#' model_fa <- "visual  =~ x1 + x2 + x3
#'              textual =~ x4 + x5 + x6
#'              speed   =~ x7 + x8 + x9
#'              pen() * visual  =~ x4 + x5 + x6 + x7 + x8 + x9
#'              pen() * textual =~ x1 + x2 + x3 + x7 + x8 + x9
#'              pen() * speed   =~ x1 + x2 + x3 + x4 + x5 + x6
#'              visual  ~~ 1 * visual
#'              textual ~~ 1 * textual
#'              speed   ~~ 1 * speed"
#'              
#' # fit with mcp under specified penalty levels and convexity levels
#' lslx_fa <- plsem(model = model_fa, 
#'                  data = lavaan::HolzingerSwineford1939,
#'                  penalty_method = "mcp", 
#'                  lambda_grid = seq(.01, .60, .01), 
#'                  delta_grid = c(1.5, 3.0, Inf))
#' 
#' # summarize fitting result under the penalty level selected by 'bic'
#' summary(lslx_fa, selector = "bic")
#' 
#' @export
## \code{plsem()} is a wrapper for \code{lslx$()$fit()}. ##
plsem <- function(model, 
                  data, 
                  penalty_method = "mcp",
                  lambda_grid = "default",
                  delta_grid = "default",
                  ...) {
  r6_lslx <- lslx$new(model = model,
                      data = data,
                      ...)
  r6_lslx$fit(penalty_method = penalty_method,
              lambda_grid = lambda_grid,
              delta_grid = delta_grid,
              ...)
  return(invisible(r6_lslx))
}


#' @export
## \code{print.lslx()} is an S3 method to summarize \code{lslx} fitting results.##
print.lslx <- function(x, ...) {
  x$print()
}



#' S3 method to summarize \code{lslx} fitting results
#' 
#' \code{summary.lslx()} is an \code{S3} interface for the \code{$summarize()} method to summarize \code{lslx} fitting results. 
#' 
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#' @param lambda A \code{numeric} to specify a chosen optimal lambda value.
#' @param delta A \code{numeric} to specify a chosen optimal lambda value.
#' @param ... Other arguments. For details, please see the \code{$summarize()} method in \code{lslx}.
#' @export
## \code{summary.lslx()} is an S3 method to summarize \code{lslx} fitting results.##
summary.lslx <- function(object,
                         selector,
                         lambda,
                         delta,
                         ...) {
  object$summarize(selector = selector,
                   lambda = lambda,
                   delta = delta)
}

#' S3 method to extract parameter estimate from \code{lslx}
#' 
#' \code{coef.lslx()} is an \code{S3} interface for the \code{$extracted_coefficient()} method to extract parameter estimate from a \code{lslx} object. 
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#' @param lambda A \code{numeric} to specify a chosen optimal lambda value.
#' @param delta A \code{numeric} to specify a chosen optimal lambda value.
#' @param ... Other arguments. For details, please see the \code{$extracted_coefficient()} method in \code{lslx}.
#' @export
## \code{coef.lslx()} is an S3 method to extract model parameters from \code{lslx}##
coef.lslx <- function(object,
                      selector,
                      lambda,
                      delta,
                      ...) {
  object$extract_coefficient(selector = selector,
                             lambda = lambda,
                             delta = delta)
}

#' S3 method to extract covariance matrix of estimates from \code{lslx}
#' 
#' \code{vcov.lslx()} is an \code{S3} interface for the \code{$extracted_coefficient_acov()} method to extract covariance matrix of parameter estimate from a \code{lslx} object.
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#' @param lambda A \code{numeric} to specify a chosen optimal lambda value.
#' @param delta A \code{numeric} to specify a chosen optimal lambda value.
#' @param ... Other arguments. For details, please see the \code{$extracted_coefficient_acov()} method in \code{lslx}.
#' @export
## \code{vcov.lslx()} is an S3 method to extract asymptotic covariance matrix of parameter estimates from \code{lslx}.##
vcov.lslx <- function(object,
                      selector,
                      lambda,
                      delta,
                      ...) {
  object$extract_coefficient_acov(selector = selector,
                                  lambda = lambda,
                                  delta = delta,
                                  ...)
}

#' S3 method to extract model-implied moments from \code{lslx}
#' 
#' \code{fitted.lslx()} is an \code{S3} interface for the \code{$extracted_implied_mean()} and \code{$extracted_implied_cov()} methods to extract model-implied moments from a \code{lslx} object.
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#' @param lambda A \code{numeric} to specify a chosen optimal lambda value.
#' @param delta A \code{numeric} to specify a chosen optimal lambda value.
#' @param ... Other arguments. For details, please see the \code{$extracted_implied_mean()} and the \code{$extracted_implied_cov()} methods in \code{lslx}.
#' @export
## \code{fitted.lslx()} is an S3 method to extract model-implied moments from \code{lslx}.##
fitted.lslx <- function(object,
                        selector,
                        lambda,
                        delta,
                        ...) {
  implied_mean <- 
    object$extract_implied_mean(selector = selector,
                                lambda = lambda,
                                delta = delta,
                                ...)
  implied_cov <- 
    object$extract_implied_cov(selector = selector,
                               lambda = lambda,
                               delta = delta,
                               ...)
  return(list(mean = implied_mean,
              cov = implied_cov))
}

#' S3 method to extract residual moments from \code{lslx}
#' 
#' \code{residuals.lslx()} is an \code{S3} interface for the \code{$extracted_residual_mean()} and \code{$extracted_residual_cov()} methods to extract residuals from a \code{lslx} object.
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#' @param lambda A \code{numeric} to specify a chosen optimal lambda value.
#' @param delta A \code{numeric} to specify a chosen optimal lambda value.
#' @param ... Other arguments. For details, please see the \code{$extracted_residual_mean()} and the \code{$extracted_residual_cov()} methods in \code{lslx}.
#' @export
## \code{residuals.lslx()} is an S3 method to extract \residual moments from \code{lslx}.##
residuals.lslx <- function(object,
                           selector,
                           lambda,
                           delta,
                           ...) {
  residual_mean <- 
    object$extract_residual_mean(selector = selector,
                                 lambda = lambda,
                                 delta = delta,
                                 ...)
  residual_cov <- 
    object$extract_residual_cov(selector = selector,
                                lambda = lambda,
                                delta = delta,
                                ...)
  return(list(mean = residual_mean,
              cov = residual_cov))
}
