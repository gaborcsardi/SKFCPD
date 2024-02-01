\name{plot_SKFCPD}
\alias{plot_SKFCPD}
\title{
Plot for SKFCPD model
}
\description{
Function to make plots on SKFCPD models after the SKFCPD model has been constructed.
}
\usage{
plot_SKFCPD(x)
}
\arguments{
  \item{x}{ an object of  class \code{SKFCPD}.}

}


\value{

Two plots: (1) plot of data with the red dashed lines mark the estimated changepoint locations, and (2) plot of the run length posterior distribution matrix. For multidimensional data, only the first dimension is plotted.
}
\references{
Li, Hanmo, Yuedong Wang, and Mengyang Gu. \emph{Sequential Kalman filter for fast online changepoint detection in longitudinal health records.} arXiv preprint arXiv:2310.18611 (2023).
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
  
  ## Detect changepoints by SKFCPD
  Online_CPD_1 = SKFCPD(design = x_mat,
                        response = y_mat,
                        train_prop = 1/3)
  
  ## visulize the results
  plot_SKFCPD(Online_CPD_1)
}
