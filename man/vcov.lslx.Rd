% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lslx-s3-interface.R
\name{vcov.lslx}
\alias{vcov.lslx}
\title{S3 method to extract covariance matrix of estimates from \code{lslx}}
\usage{
\method{vcov}{lslx}(object, selector, lambda, delta, ...)
}
\arguments{
\item{object}{A fitted \code{lslx} object.}

\item{selector}{A \code{character} to specify a selector for determining an optimal penalty level.
Its value can be any one in \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
or their robust counterparts \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"} if raw data is available.}

\item{lambda}{A \code{numeric} to specific a chosen optimal penalty level. 
If the specified \code{lambda} is not in \code{lambda_grid}, a nearest legitimate value will be used.}

\item{delta}{A \code{numeric} to specific a chosen optimal convexity level.
If the specified \code{delta} is not in \code{delta_grid}, a nearest legitimate value will be used.}

\item{...}{Other arguments. For details, please see the \code{$extracted_coefficient_acov()} method in \code{lslx}.}
}
\description{
\code{vcov.lslx()} is an \code{S3} interface for extracting covariance matrix of parameter estimate from a \code{lslx} object.
}
