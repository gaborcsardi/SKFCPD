\name{SKFCPD-package}
\alias{SKFCPD-package}
\docType{package}
\title{
  Dynamic Linear Model for Online Changepoint Detection
}
\description{
  The `SKFCPD` package provides estimation of changepoint locations using the Dynamic Linear Model (DLM) within the Bayesian Online Changepoint Detection (BOCPD) framework. The efficient computation is achieved through implementation of the Sequential Kalman filter. The range parameter and noise-to-signal ratio are estimated from training samples via a Gaussian process model. This package is capable of handling multidimensional data with temporal correlations and random missing patterns.
}
\details{
  The DESCRIPTION file:
\packageDESCRIPTION{SKFCPD}
\packageIndices{SKFCPD}
Implements a fast online changepoint detection algorithm using dynamic linear model based on Sequential Kalman filter. It's for temporally correlated data and accepts multi-dimensional datasets with missing values.
}
\author{
\packageAuthor{SKFCPD}

Maintainer: \packageMaintainer{SKFCPD}
}
\references{

Li, Hanmo, Yuedong Wang, and Mengyang Gu. \emph{Sequential Kalman filter for fast online changepoint detection in longitudinal health records.} arXiv preprint arXiv:2310.18611 (2023).

Fearnhead, Paul, and Zhen Liu. \emph{On-line inference for multiple changepoint problems.} Journal of the Royal Statistical Society Series B: Statistical Methodology 69, no. 4 (2007): 589-605.

Adams, Ryan Prescott, and David JC MacKay. \emph{Bayesian online changepoint detection.} arXiv preprint arXiv:0710.3742 (2007).

Hartikainen, Jouni, and Simo Sarkka. \emph{Kalman filtering and smoothing solutions to temporal Gaussian process regression models.} In 2010 IEEE international workshop on machine learning for signal processing, pp. 379-384. IEEE, 2010.


Gu, Mengyang, and Yanxun Xu. \emph{Fast nonseparable Gaussian stochastic process with application to methylation level interpolation.} Journal of Computational and Graphical Statistics 29, no. 2 (2020): 250-260.

Gu, Mengyang, and Weining Shen. \emph{Generalized probabilistic principal component analysis of correlated data.} The Journal of Machine Learning Research 21, no. 1 (2020): 428-468.

Gu, Mengyang, Xiaojing Wang, and James O. Berger. \emph{Robust Gaussian stochastic process emulation.} The Annals of Statistics 46, no. 6A (2018): 3038-3066.
}

\seealso{
 \code{\link[SKFCPD]{SKFCPD}}
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
}
