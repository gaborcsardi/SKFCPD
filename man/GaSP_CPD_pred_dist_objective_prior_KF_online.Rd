\name{GaSP_CPD_pred_dist_objective_prior_KF_online}
\alias{GaSP_CPD_pred_dist_objective_prior_KF_online}

\title{ Computing the predictive distribution in the online fashion
}
\description{
This function computs the predictive distribution of the run length in the online fashion.
}
\usage{
  GaSP_CPD_pred_dist_objective_prior_KF_online(KF_params, prev_L_params, cur_point,
  d, gamma, model_type, mu, sigma_2, eta, kernel_type, G_W_W0_V_ini, G_W_W0_V)
}

\arguments{
    \item{KF_params}{A list of current Kalman filter parameters.}
    \item{prev_L_params}{A list of previous Kalman filter parameters.}
    \item{cur_point}{A value of current observation.}
    \item{d}{A value of the distance between the sorted input.}
    \item{gamma}{A numeric variable of the range parameter for the covariance matrix. The default value of gamma is 1.}
    \item{model_type}{A numeric variable that can take values of 0, 1 and 2. Model_type=0 stands for a GP model with unknown mean and known variance. Model_type=1 stands for a GP model with known mean and unknown variance. Model_type=2 stands for a GP model with unknown mean and unknown variance. The default value of model_type is 2.}
    \item{mu}{A vector of the mean parameter at each coordinate. Ignored when model_type = 0 or 2.}
    \item{sigma_2}{A vector of the variance parameter at each coordinate.}
    \item{eta}{A vector of the noise-to-signal ratio at each coordinate}
    \item{kernel_type}{A character specifying the type of kernels of the input. \code{matern_5_2} are \code{Matern correlation} with roughness parameter 5/2. \code{exp} is power exponential correlation with roughness parameter alpha=2. The default choice is \code{matern_5_2}.}
    \item{G_W_W0_V_ini}{A list of the initial coefficient and conditional matrics for Gaussian Process(GP) model. It's the output from the function \code{Construct_G_W_W0_V} }
    \item{G_W_W0_V}{A list of the coefficient and conditional matrics for Gaussian Process(GP) model. It's the output from the function \code{Construct_G_W_W0_V} }
}
\value{
\code{GaSP_CPD_pred_dist_objective_prior_KF_online} returns a list that contains 3 items: (1) the current Kalman filter parameters; (2) the previous Kalman filter parameters and (3) the vector of the logrithm for the current predictive distribution of different run lengths.
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