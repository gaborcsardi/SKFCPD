\name{KF_param_update_for_profile_like}
\alias{KF_param_update_for_profile_like}

\title{ Updating Kalman filter parameters
}
\description{
Updating the Kalman filter parameters for Gaussian Process model with Matern 2.5 or power exponential kernels with different observation sequences.
}
\usage{
  KF_param_update_for_profile_like(cur_input, cur_num_obs,
  prev_param, d, gamma, eta, kernel_type, G_W_W0_V)
}

\arguments{
    \item{cur_input}{A value of current observation.}
    \item{cur_num_obs}{A value of index for the current observation.}
    \item{prev_param}{A list of previous Kalman filter parameters.}
    \item{d}{A value of the distance between the sorted input.}
    \item{gamma}{A value of the range parameter for the covariance matrix.}
    \item{eta}{The noise-to-signal ratio.}
    \item{kernel_type}{A character specifying the type of kernels of the input. \code{matern_5_2} are \code{Matern correlation} with roughness parameter 5/2. \code{exp} is power exponential correlation with roughness parameter alpha=2.}
    \item{G_W_W0_V}{A list of the coefficient and conditional matrics for Gaussian Process(GP) model. It's the output from the function \code{Construct_G_W_W0_V} }
}
\value{
\code{KF_param_update_for_profile_like} returns a list of updated kalman filter parameters with different observation sequences.
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