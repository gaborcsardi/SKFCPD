\name{get_predictive_dist_direct_objective_prior}
\alias{get_predictive_dist_direct_objective_prior}

\title{ Updating the predictive distribution
}
\description{
Updating the predictive distribution of the run length under the objective prior directly.
}
\usage{
  get_predictive_dist_direct_objective_prior(cur_input_seq, d, gamma, mu, sigma_2, eta)
}

\arguments{
    \item{cur_input_seq}{A vector of sequence of observations.}
    \item{d}{A value of the distance between the sorted input.}
    \item{gamma}{A numeric variable of the range parameter for the covariance matrix. The default value of gamma is 1.}
    \item{eta}{A vector of the noise-to-signal ratio at each coordinate}
    \item{mu}{A vector of the mean parameter at each coordinate. Ignored when model_type = 0 or 2.}
    \item{sigma_2}{A vector of the variance parameter at each coordinate.}
}
\value{
\code{get_predictive_dist_direct_objective_prior} returns the log likelihood of observations that follows Gaussian Process with Exponential kernel.
}

\references{
Williams, C. K., & Rasmussen, C. E. (2006). \emph{Gaussian processes for machine learning (Vol. 2, No. 3, p. 4).} \emph{Cambridge, MA: MIT press.}
}
\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\keyword{internal}