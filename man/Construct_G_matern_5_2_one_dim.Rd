\name{Construct_G_matern_5_2_one_dim}
\alias{Construct_G_matern_5_2_one_dim}
\title{ 
The coefficient matrix in the dynamic linear model when kernel is the Matern covariance with roughness parameter 2.5. 
}
\description{
The coefficient matrix in the dynamic linear model when kernel is the Matern covariance with roughness parameter 2.5. 
}
\usage{
Construct_G_matern_5_2_one_dim(delta_x,lambda)
}
\arguments{
  \item{delta_x}{A value of the distance between the sorted input.
}
  \item{lambda}{the transformed range parameter. }
}


\value{
GG matrix. 
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
