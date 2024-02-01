\name{Construct_W_exp_fastGP}
\alias{Construct_W_exp_fastGP}
\title{ 
The conditional covariance matrix of the state in the dynamic linear model when kernel is the exponential covariance
}
\description{
The conditional covariance matrix of the state in the dynamic linear model when kernel is the exponential covariance. 
}
\usage{
Construct_W_exp_fastGP(delta_x,lambda,W0)
}
\arguments{

  \item{delta_x}{the distance between the sorted input.
}
  \item{lambda}{the transformed range parameter. }
    \item{W0}{the covariance matrix of the stationary distribution of the state. }

}


\value{
W matrix. 
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
