\name{CPD_DLM}
\alias{CPD_DLM}

\title{ Setting up the CPD_DLM model
}
\description{
Implementing the robust GaSP model for estimating the changepoint locations. The range parameter and noise-to-signal ratio are estimated from the training samples by a Gaussian process model.
}
\usage{
  CPD_DLM(design, response, gamma,model_type, mu, sigma_2, eta,
         kernel_type, stop_at_first_cp, hazard_vec,
         truncate_at_prev_cp)
}

\arguments{
    \item{design}{A matrix with dimension n x p. The design of the experiment. }
    \item{response}{A matrix with dimension n x q. The observations. }
    \item{gamma}{A numeric variable of the range parameter for the covariance matrix. The default value of gamma is 1.}
    \item{model_type}{A numeric variable that can take values of 0, 1 and 2. Model_type=0 stands for a GP model with unknown mean and known variance. Model_type=1 stands for a GP model with known mean and unknown variance. Model_type=2 stands for a GP model with unknown mean and unknown variance. The default value of model_type is 2.}
    \item{mu}{A vector of the mean parameter at each coordinate. Ignored when model_type = 0 or 2.}
    \item{sigma_2}{A vector of the variance parameter at each coordinate.}
    \item{eta}{A vector of the noise-to-signal ratio at each coordinate}
    \item{kernel_type}{A character specifying the type of kernels of the input. \code{matern_5_2} are \code{Matern correlation} with roughness parameter 5/2. \code{exp} is power exponential correlation with roughness parameter alpha=2. The default choice is \code{matern_5_2}.}
    \item{stop_at_first_cp}{A numeric variable that decides if the SKFCPD method stops when it detects the first changepoint. The default value of stop_at_first_cp is FALSE.}
    \item{hazard_vec}{The hazard vector in the SKFCPD method. 1/vector is the prior probability that a changepoint occur at a vector of time points.}
    \item{truncate_at_prev_cp}{If TRUE, truncate the run length at the most recently detected changepoint. The default value of truncate_at_prev_cp is FALSE.}
}
\value{
\code{SKFCPD} returns a S4 object of class \code{SKFCPD} (see \code{SKFCPD-class}).
}

\references{
Li, Hanmo, Yuedong Wang, and Mengyang Gu. \emph{Sequential Kalman filter for fast online changepoint detection in longitudinal health records.} arXiv preprint arXiv:2310.18611 (2023).

Fearnhead, Paul, and Zhen Liu. \emph{On-line inference for multiple changepoint problems.} Journal of the Royal Statistical Society Series B: Statistical Methodology 69, no. 4 (2007): 589-605.

Adams, Ryan Prescott, and David JC MacKay. \emph{Bayesian online changepoint detection.} arXiv preprint arXiv:0710.3742 (2007).

Hartikainen, Jouni, and Simo Sarkka. \emph{Kalman filtering and smoothing solutions to temporal Gaussian process regression models.} In 2010 IEEE international workshop on machine learning for signal processing, pp. 379-384. IEEE, 2010.
}
\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}

\keyword{internal}