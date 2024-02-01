\name{get_LY_online}
\alias{get_LY_online}

\title{ Updating Kalman filter parameters
}
\description{
Updating the Kalman filter parameters for Gaussian Process model with Matern 2.5 or power exponential kernels.
}
\usage{
  get_LY_online(cur_input, prev_param, eta, G_W_W0_V)
}

\arguments{
    \item{cur_input}{A value of current observation.}
    \item{prev_param}{A list of previous Kalman filter parameters.}
    \item{eta}{The noise-to-signal ratio.}
    \item{G_W_W0_V}{A list of the coefficient and conditional matrics for Gaussian Process(GP) model. It's the output from the function \code{Construct_G_W_W0_V} }
}
\value{
\code{get_LY_online} returns a list of updated kalman filter parameters.
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