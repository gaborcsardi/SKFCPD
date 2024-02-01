\name{Get_Q_K}
\alias{Get_Q_K}
\title{ 
 matrices and vectors for the inverse covariance in the predictive distribution
}
\description{
This function computes the required values for the inverse covariance matrix. 
}
\usage{
Get_Q_K(GG,W,C0,VV)
}
\arguments{
    \item{GG}{a list of matrices defined in the dynamic linear model.}
    \item{W}{a list of matrices defined in the dynamic linear model.}
    \item{C0}{a matrix defined in the dynamic linear model.}
    \item{VV}{a numerical value for the nugget.}
}


\value{
A list of 2 items for Q and K. 
}

\references{

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models}.  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2019), \emph{fast nonseparable gaussian stochastic process with application to methylation level interpolation}.  \emph{Journal of Computational and Graphical Statistics}, In Press, arXiv:1711.11501.

Campagnoli P, Petris G, Petrone S. (2009), \emph{Dynamic linear model with R}. Springer-Verlag New York.

}

\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\keyword{internal}
