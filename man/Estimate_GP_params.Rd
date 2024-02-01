\name{Estimate_GP_params}
\alias{Estimate_GP_params}

\title{Estimate parameters from fast computation of GaSP model
}
\description{
Getting the estimated parameters from fast computation of the Gaussian stochastic process (GaSP) model with the Matern kernel function with a noise.
}
\usage{
  Estimate_GP_params(input, output, kernel_type='matern_5_2')
}

\arguments{
    \item{input}{a vector with dimension num_obs x 1 for the sorted input locations.}
    \item{output}{a  vector with dimension n x 1 for the observations at the sorted input locations.}
  \item{kernel_type}{a \code{character} to specify the type of kernel to use. The current version supports kernel_type to be "matern_5_2" or "exp", meaning that the matern kernel with roughness parameter being 2.5 or 0.5 (power exponent kernel), respectively. 
}
}
\value{

\code{Estimate_GP_params} returns an S4 object of class \code{Estimated_GP_params} with estimated parameters including

  \item{beta }{the inverse range parameter, i.e. beta=1/gamma}
  \item{eta }{the noise-to-signal ratio}
  \item{sigma_2 }{the variance parameter}
}

\references{
Hartikainen, Jouni, and Simo Sarkka. \emph{Kalman filtering and smoothing solutions to temporal Gaussian process regression models.} In 2010 IEEE international workshop on machine learning for signal processing, pp. 379-384. IEEE, 2010.


Gu, Mengyang, and Yanxun Xu. \emph{Fast nonseparable Gaussian stochastic process with application to methylation level interpolation.} Journal of Computational and Graphical Statistics 29, no. 2 (2020): 250-260.

Gu, Mengyang, and Weining Shen. \emph{Generalized probabilistic principal component analysis of correlated data.} The Journal of Machine Learning Research 21, no. 1 (2020): 428-468.

Gu, Mengyang, Xiaojing Wang, and James O. Berger. \emph{Robust Gaussian stochastic process emulation.} The Annals of Statistics 46, no. 6A (2018): 3038-3066.
}
\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\examples{

  library(SKFCPD)

  #------------------------------------------------------------------------------
  # simple example with noise
  #------------------------------------------------------------------------------
  
  y_R<-function(x){
    cos(2*pi*x)
  }
  ###let's test for 100 observations
  set.seed(1)
  num_obs=100
  input=runif(num_obs)
  output=y_R(input)+rnorm(num_obs,mean=0,sd=1)
  
  ## run Estimate_GP_params to get estimated parameters
  params_est = Estimate_GP_params(input, output)
  print(params_est@beta) ## inverse of range parameter
  print(params_est@eta) ## noise-to-signal ratio
  print(params_est@sigma_2) ## variance
}
