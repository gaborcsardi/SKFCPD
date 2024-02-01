\name{Construct_W0_exp_one_dim}
\alias{Construct_W0_exp_one_dim}
\title{ 
covariance of the stationary distribution of the state when kernel is the exponential covariance.
}
\description{
This function computes the covariance of the stationary distribution of the state when kernel is the exponential covariance.
}
\usage{
Construct_W0_exp_one_dim(lambda)
}
\arguments{

  \item{lambda}{the transformed range parameter. }
}


\value{
W0 matrix. 
}

\references{

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models}.  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2019), \emph{fast nonseparable Gaussian stochastic process with application to methylation level interpolation}.  \emph{Journal of Computational and Graphical Statistics}, In Press, arXiv:1711.11501.

Campagnoli P, Petris G, Petrone S. (2009), \emph{Dynamic linear model with R}. Springer-Verlag New York.

}

\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\keyword{internal}
