\name{compute_log_lik}
\alias{compute_log_lik}
\title{ 
Natural logarithm of profile likelihood by the fast computing algorithm
}
\description{
This function computes the natural logarithm of the profile likelihood for the range and nugget parameter after plugging the closed form maximum likelihood estimator for the variance parameter. 
}
\usage{
compute_log_lik(param, design, response, kernel_type)
}
\arguments{
  \item{param}{a vector of parameters. The first parameter is the natural logarithm of the inverse range parameter in the kernel function. If the data contain noise, the second parameter is the logarithm of the nugget-variance ratio parameter.}
  \item{design}{A matrix with dimension n x p. The design of the experiment. }
  \item{response}{A matrix with dimension n x q. The observations. }
    \item{kernel_type}{A character specifying the type of kernels of the input. \code{matern_5_2} are \code{Matern correlation} with roughness parameter 5/2. \code{exp} is power exponential correlation with roughness parameter alpha=2. The default choice is \code{matern_5_2}.}
}

\value{
The numerical value of natural logarithm of the profile likelihood.
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

\keyword{internal}