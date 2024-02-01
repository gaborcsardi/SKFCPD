\name{get_predictive_dist_KF_objective_prior}
\alias{get_predictive_dist_KF_objective_prior}

\title{ Updating the predictive distribution
}
\description{
Updating the predictive distribution of the run length under the objective prior.
}
\usage{
  get_predictive_dist_KF_objective_prior(cur_input, cur_num_obs,
  params, prev_L, d, gamma, model_type, mu, sigma_2, eta, kernel_type)
}

\arguments{
    \item{cur_input}{A value of current observation.}
    \item{cur_num_obs}{A value of index for the current observation.}
    \item{params}{A list of current Kalman filter parameters.}
    \item{prev_L}{A list of previous Kalman filter parameters.}
    \item{d}{A value of the distance between the sorted input.}
    \item{gamma}{A numeric variable of the range parameter for the covariance matrix.}
    \item{model_type}{A numeric variable that can take values of 0, 1 and 2. Model_type=0 stands for a GP model with unknown mean and known variance. Model_type=1 stands for a GP model with known mean and unknown variance. Model_type=2 stands for a GP model with unknown mean and unknown variance.}
    \item{mu}{A vector of the mean parameter at each coordinate. Ignored when model_type = 0 or 2.}
    \item{sigma_2}{A vector of the variance parameter at each coordinate.}
    \item{eta}{A vector of the noise-to-signal ratio at each coordinate}
    \item{kernel_type}{A character specifying the type of kernels of the input. \code{matern_5_2} are \code{Matern correlation} with roughness parameter 5/2. \code{exp} is power exponential correlation with roughness parameter alpha=2.}
}
\value{
\code{get_predictive_dist_KF_objective_prior} returns a list of updated predictive distribution of the run length under the objective prior.
}

\references{
Fearnhead, P., & Liu, Z. (2007). \emph{On-line inference for multiple changepoint problem.} \emph{ Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 69(4), 589-605.

Adams, R. P., & MacKay, D. J. (2007). \emph{Bayesian online changepoint detection.} \emph{arXiv preprint} arXiv:0710.3742.

  Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models},  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.
}
\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\keyword{internal}