\name{SKFCPD-class}
\Rdversion{1.1}
\docType{class}
\alias{SKFCPD-class}

\title{Class \code{"SKFCPD"}}
\description{
S4 class for SKFCPD where the range parameter and noise-to-signal parameters are estimated from the training samples.
}
\section{Objects from the Class}{
Objects of this class are created and initialized with the function \code{\link{SKFCPD}} that computes the calculations needed for setting up the analysis.
}
\section{Slots}{
  \describe{
    \item{\code{design}:}{Object of class \code{"matrix"} with dimension n x p. The design of the experiment. }
    \item{\code{response}:}{Object of class \code{"matrix"} with dimension n x q. The observations. }
    \item{\code{test_start}:}{Object of class \code{"numeric"}. The starting index of test period.}
    \item{\code{kernel_type}:}{Object of class \code{"character"} to specify the type of kernel to use. }
    \item{\code{gamma}:}{Object of class \code{"vector"} with dimension q x 1. The range parameters. }
    \item{\code{eta}:}{Object of class \code{"vector"} with dimension q x 1. The noise-to-signal ratio. }
    \item{\code{sigma_2}:}{Object of class \code{"vector"} with dimension q x 1. The variance parameters. }
    \item{\code{hazard_vec}:}{Object of class \code{"numeric"}. The n x 1 hazard vector in the FastCPD method. }
    \item{\code{KF_params_list}:}{Object of class \code{"list"}. The list of Kalman filter parameters from the previous run of the algorithm. }
    \item{\code{prev_L_params_list}:}{Object of class \code{"list"}. The list of parameters for calculating the quadratic form of the inverse covariance matrix from the previous run of the algorithm. }
    \item{\code{run_length_posterior_mat}:}{Object of class \code{"matrix"} with dimension n x n. The posterior distribution of the run length. }
    \item{\code{run_length_joint_mat}:}{Object of class \code{"matrix"} with dimension n x n. The joint distribution of the run length and the observations. }
    \item{\code{log_pred_dist_mat}:}{Object of class \code{"matrix"} with dimension n x n. The logrithm of the predictive distribution of observations. }
    \item{\code{cp}:}{Object of class \code{"vector"} with length m. The location of estimated changepoints.}
  }
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

\seealso{ \code{\link{SKFCPD}} for more details about how to create a \code{SKFCPD} object.
}

\keyword{classes}
