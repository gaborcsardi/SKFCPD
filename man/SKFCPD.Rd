\name{SKFCPD}
\alias{SKFCPD}

\title{ Getting the results of the SKFCPD model
}
\description{
Estimating changepoint locations using the Dynamic Linear Model (DLM) within the Bayesian Online Changepoint Detection (BOCPD) framework. The efficient computation is achieved through implementation of the Kalman filter. The range parameter and noise-to-signal ratio are estimated from training samples via a Gaussian process model. This function is capable of handling multidimensional data with temporal correlations and random missing patterns.
}
\usage{
  SKFCPD(design = NULL, response = NULL, FCPD = NULL, 
  init_params = list(gamma = 1, sigma_2 = 1, eta = 1), 
  train_prop = NULL, kernel_type = "matern_5_2", 
  hazard_vec=100, print_info = TRUE, truncate_at_prev_cp = FALSE)
}

\arguments{
    \item{design}{A vector with the length of n. The design of the experiment. }
    \item{response}{A matrix with dimension n x q. The observations. }
    \item{FCPD}{An object of the class \code{SKFCPD} computed in the previous run of the algorithm.}
    \item{init_params}{A list with estimated range parameter \code{gamma}, noise-to-signal parameter \code{eta} and variance parameter \code{sigma_2}. The default values are \code{gamma}=1, \code{eta}=1, and \code{sigma_2}=1.}
    \item{train_prop}{A numerical value between 0 and 1. The propotation of training samples for parameter estimation. When \code{train_prop}=NULL, we skip the training process and specify the parameter values in the argument \code{init_params}.}
    \item{kernel_type}{A character specifying the type of kernels of the input. \code{matern_5_2} are Matern correlation with roughness parameter 5/2. \code{exp} is power exponential correlation with roughness parameter alpha=2. The default choice is \code{matern_5_2}.}
    \item{hazard_vec}{Either a constant or a vector with the length of n. The hazard vector in the SKFCPD method. hazard_vec = 1/\code{hazard_const} is the prior probability that a changepoint occur at any time points. The default value of hazard_vec is 100. }
    \item{print_info}{This setting prints out updates on the progress of the algorithm if set to TRUE.}
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

\examples{
  library(SKFCPD)
  
  #------------------------------------------------------------------------------
  # Example: fast online changepoint detection with DEPENDENT data.
  # 
  # Data generation: Data follows a multidimensional Gaussian process with Matern 2.5 kernel.
  #------------------------------------------------------------------------------
  # Data Generation
  set.seed(1)
  
  n_obs = 150
  n_dim = 2
  seg_len = c(70, 30, 20,30)
  mean_each_seg = c(0,1,-1,0)
  
  x_mat=matrix(1:n_obs)
  y_mat=matrix(NA, nrow=n_obs, ncol=n_dim)
  
  gamma = rep(5, n_dim) # range parameter of the covariance matrix
  
  # compute the matern 2.5 kernel
  construct_cor_matrix = function(input, gamma){
    n = length(input)
    R0=abs(outer(input,(input),'-'))
    matrix_one = matrix(1, n, n)
    const = sqrt(5) * R0 / gamma
    Sigma = (matrix_one + const + const^2/3) * (exp(-const))
    return(Sigma)
  }
  
  for(j in 1:n_dim){
    y_each_dim = c()
    for(i in 1:length(seg_len)){
      nobs_per_seg = seg_len[i]
      Sigma = construct_cor_matrix(1:nobs_per_seg, gamma[j])
      L=t(chol(Sigma))
      theta=rep(mean_each_seg[i],nobs_per_seg)+L\%*\%rnorm(nobs_per_seg)
      y_each_dim = c(y_each_dim, theta+0.1*rnorm(nobs_per_seg))
    }
    y_mat[,j] = y_each_dim
  }
  
  plot(y_mat[,1]~x_mat)
  
  ## Detect changepoints by SKFCPD
  Online_CPD_1 = SKFCPD(design = x_mat,
                        response = y_mat,
                        train_prop = 1/3)
  
  ## visulize the results
  plot_SKFCPD(Online_CPD_1)
  
  #------------------------------------------------------------------------------
  # Example: iterative fast online changepoint detection
  #------------------------------------------------------------------------------
  
  # train/test split
  train_index = 1:50
  test_index = 51:150
  
  x_train = x_mat[train_index, , drop=FALSE]
  y_train = y_mat[train_index, , drop=FALSE]
  x_test = x_mat[test_index, , drop=FALSE]
  y_test = y_mat[test_index, , drop=FALSE]
  
  # estimate the parameters
  fast_model = Estimate_GP_params(x_train, y_train) 
  
  # detect changepoints sequentially
  Online_CPD_2 = NULL
  init_params = list(gamma = 1/fast_model@beta,
                     eta = fast_model@eta,
                     sigma_2 = fast_model@sigma_2)
  for(i in 1:length(test_index)){
    if(sum(is.na(y_test[i,]))==0){
      Online_CPD_2 = SKFCPD(design = x_test[i, , drop=FALSE],
                            response = y_test[i, , drop=FALSE],
                            FCPD = Online_CPD_2,
                            init_params = init_params)
    }
  }
}
