\name{Estimated_GP_params-class}
\docType{class}
\alias{Estimated_GP_params-class}
\title{Estimated GaSP parameters class}
\description{S4 class for fast parameter estimation of the Gaussian stochastic process (GaSP) model with the Matern kernel function with or without a noise.}

\section{Objects from the Class}{
Objects of this class are created with the function \code{\link{Estimate_GP_params}} that computes the calculations needed for setting up the estimation and prediction.}

\section{Slots}{
	 \describe{
    \item{\code{beta}:}{object of class \code{numeric} for the inverse of the range parameter, i.e. beta = 1/gamma.}
    \item{\code{eta}:}{object of class \code{numeric} for the estimated noise-to-signal parameter.}
    \item{\code{sigma_2}:}{object of class \code{numeric} for the estimated variance parameter.}
    }
}

\references{

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models},  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2017), \emph{Nonseparable Gaussian stochastic process: a unified
view and computational strategy}, arXiv:1711.11501.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian Stochastic Process Emulation}, \emph{Annals of Statistics}, \bold{46}, 3038-3066.

}

\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\keyword{classes}
\keyword{internal}
