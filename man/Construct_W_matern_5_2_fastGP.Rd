\name{Construct_W_matern_5_2_fastGP}
\alias{Construct_W_matern_5_2_fastGP}
\title{ 
The conditional covariance matrix for  matern covariance with roughness parameter 2.5
}
\description{
The conditional covariance matrix of the state in the dynamic linear model when kernel is the matern covariance with roughness parameter 2.5.
}
\usage{
Construct_W_matern_5_2_fastGP(delta_x,lambda, W0)
}
\arguments{
  \item{delta_x}{a vector of the distance between the sorted input.
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
