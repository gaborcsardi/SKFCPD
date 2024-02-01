\name{Construct_G_W_W0_V}
\alias{Construct_G_W_W0_V}
\title{ 
Generating coefficient and conditional matrics}
\description{
Generating coefficient and conditional matrics for Gaussian Process(GP) model with Matern 2.5 or power exponential kernels.
}
\usage{
Construct_G_W_W0_V(d, gamma, eta, kernel_type, is_initial)
}
\arguments{
  \item{d}{A value of the distance between the sorted input.}
  \item{gamma}{A value of the range parameter for the covariance matrix.}
  
  \item{eta}{The noise-to-signal ratio. }
  
  \item{kernel_type}{A character specifying the type of kernels of the input. \code{matern_5_2} are \code{Matern correlation} with roughness parameter 5/2. \code{exp} is power exponential correlation with roughness parameter alpha=2. The default choice is \code{matern_5_2}.}
  
  \item{is_initial}{A bolean variable. is_initial=TRUE means the matrics generated is for the inital state.}
}


\value{
A list of GG, W, W0 and VV matrix. 
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
