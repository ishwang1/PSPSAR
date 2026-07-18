\name{SAR}
\alias{SAR}
\title{Estimate a Constant Coefficient Spatial Autoregressive Model}
\usage{
SAR(Y, X, W, iv_lag = 1, estimator = c("2SLS", "GMM"),
  ci_level = 0.95, hac_method = c("Newey-West", "Andrews"),
  preprocess_W = FALSE)
}
\arguments{
  \item{Y}{A numeric \code{T} by \code{N} outcome matrix.}
  \item{X}{An \code{N} by \code{p} numeric time-invariant covariate matrix if not time-varying, or a \code{T} by \code{N} by \code{p} numeric time-varying array. An intercept is added internally.}
  \item{W}{An \code{N} by \code{N} spatial weight matrix.}
  \item{iv_lag}{A positive integer giving the maximum spatial lag of covariates included as instruments.}
  \item{estimator}{\code{"2SLS"} or \code{"GMM"}.}
  \item{ci_level}{The confidence level for confidence intervals.}
  \item{hac_method}{\code{"Newey-West"} for Newey-West HAC or \code{"Andrews"} for Andrews HAC.}
  \item{preprocess_W}{If \code{FALSE}, apply built-in spatial weight preprocessing when needed. If \code{TRUE}, use \code{W} as supplied.}
}
\value{A list with \code{beta}, \code{rho}, \code{residuals}, and \code{diagnostics}. Coefficient tables contain \code{term}, \code{estimate}, \code{se}, \code{lower}, and \code{upper}; diagnostics include estimation, first-stage, HAC, and spatial weight details.}
\details{The model has one constant spatial autoregressive coefficient. It instruments the endogenous spatial lag using covariates and their spatial lags, and estimates HAC standard errors over time.}
\description{Fits a constant coefficient spatial autoregressive benchmark using two-stage least squares or feasible GMM.}
