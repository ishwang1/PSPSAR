\name{PSPSAR}
\alias{PSPSAR}
\title{Estimate a Profiled State-Varying Spatial Autoregressive Model}
\usage{
PSPSAR(Y, X, W, z, grid = sort(unique(z)), iv_lag = 1,
  estimator = c("2SLS", "GMM"), inference = c("debiased", "undersmooth"),
  bandwidth = NULL, pilot_bandwidth = NULL,
  kernel = c("gaussian", "epanechnikov", "triangular", "uniform", "quartic"),
  ci_level = 0.95, hac_method = c("Newey-West", "Andrews"),
  preprocess_W = FALSE)
}
\arguments{
  \item{Y}{A numeric \code{T} by \code{N} outcome matrix; rows are periods and columns are spatial units.}
  \item{X}{An \code{N} by \code{p} numeric time-invariant covariate matrix, or a \code{T} by \code{N} by \code{p} numeric time-varying array. An intercept is added internally.}
  \item{W}{An \code{N} by \code{N} spatial weight matrix.}
  \item{z}{A numeric state vector of length \code{T}.}
  \item{grid}{Numeric state values at which to evaluate the spillover curve.}
  \item{iv_lag}{A positive integer giving the maximum spatial lag of covariates included as instruments. Spatial lags are generated only for covariates with cross-sectional variation in at least one period.}
  \item{estimator}{\code{"2SLS"} for profiled two-stage least squares or \code{"GMM"} for feasible profiled GMM.}
  \item{inference}{\code{"debiased"} for local quadratic bias correction or \code{"undersmooth"} for uncorrected inference.}
  \item{bandwidth}{An optional positive main smoothing bandwidth.}
  \item{pilot_bandwidth}{An optional positive pilot bandwidth for debiased inference.}
  \item{kernel}{The kernel used for local smoothing.}
  \item{ci_level}{The confidence level for confidence intervals.}
  \item{hac_method}{\code{"Newey-West"} for Newey-West HAC or \code{"Andrews"} for Andrews HAC.}
  \item{preprocess_W}{If \code{FALSE}, apply built-in diagonal removal, isolated unit removal, and row normalization when needed. If \code{TRUE}, use \code{W} as supplied.}
}
\value{A list with \code{beta}, \code{rho}, \code{residuals}, \code{hac_bandwidths}, and \code{diagnostics}. \code{beta} and \code{rho} contain raw estimates, debiased estimates in \code{debiased}, standard errors, confidence intervals, and p-values in \code{p}; \code{rho} is evaluated at \code{grid}. The diagnostics identify covariates included in and excluded from the spatial-lag instrument set.}
\details{The model has a state-varying spatial coefficient and global covariate coefficients. It accepts time-invariant or time-varying covariates and uses spatial lags of covariates as instruments. A covariate that is common to all spatial units within every period remains in the outcome equation but does not receive a redundant spatial-lag instrument.}
\references{Chen, Haiqiang, Yingxing Li, Xiaojun Song, and Yishu Wang (2026). \emph{Profiled Semiparametric GMM for a State-Varying SAR Coefficient}. Working paper.\cr\cr Lan, Zhiqiang, Xiangfu Luo, Yishu Wang, and Guoyao Wu (2026). \emph{Air Quality and Spatial Spillovers in Electric Vehicle Charging: A State-Varying Spatial Autoregressive Model}. Working paper.}
\description{Fits a profiled state-varying spatial autoregressive model with 2SLS or feasible GMM estimation and HAC inference.}
