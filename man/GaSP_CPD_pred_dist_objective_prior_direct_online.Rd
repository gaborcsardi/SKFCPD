\name{GaSP_CPD_pred_dist_objective_prior_direct_online}
\alias{GaSP_CPD_pred_dist_objective_prior_direct_online}

\title{ Computing the predictive distribution directly in the online fashion
}
\description{
This function computs directly the predictive distribution of the run length in the online fashion. The direct computation includes the inversion of covariance matrix, which is of computational complexity $O(n^3)$, with $n$ being the number of observations.
}
\usage{
  GaSP_CPD_pred_dist_objective_prior_direct_online(cur_seq, d, gamma, eta, mu, sigma_2)
}

\arguments{
    \item{cur_seq}{A vector of sequence of observations.}
    \item{d}{A value of the distance between the sorted input.}
    \item{gamma}{A numeric variable of the range parameter for the covariance matrix. The default value of gamma is 1.}
    \item{eta}{A vector of the noise-to-signal ratio at each coordinate}
    \item{mu}{A vector of the mean parameter at each coordinate. Ignored when model_type = 0 or 2.}
    \item{sigma_2}{A vector of the variance parameter at each coordinate.}
}
\value{
\code{GaSP_CPD_pred_dist_objective_prior_direct_online} returns the log likelihood of observations that follows Gaussian Process with Exponential kernel.
}

\references{
Williams, C. K., & Rasmussen, C. E. (2006). \emph{Gaussian processes for machine learning (Vol. 2, No. 3, p. 4).} \emph{Cambridge, MA: MIT press.}
}
\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\keyword{internal}