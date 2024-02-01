setClass("Estimated_GP_params", 		
         representation( 
           beta = "numeric",          ## inverse of the range parameter, i.e. beta = 1/gamma
           eta="numeric",            ## noise-to-variance ratio
           sigma_2="numeric"           ## variance
         ),
)



setClass("SKFCPD", 		
         representation( 
           ## data
           design = "matrix", ## location of the sort input
           response = "matrix", ## the observations, size nxp
           test_start = "numeric",
           ## parameter
           kernel_type = "character", ## information of the kernel
           gamma = "vector", ## the range parameter
           eta = "vector", ## the nugget parameter
           sigma_2 = "vector", ## the nugget parameter
           hazard_vec = "numeric", ## the hazard vector
           ## results
           KF_params_list = "list",
           prev_L_params_list = "list",
           run_length_posterior_mat = "matrix",
           run_length_joint_mat = "matrix",
           log_pred_dist_mat = "matrix",
           cp = "vector"
         ),
)


if(!isGeneric("plot")) {
  setGeneric(name = "plot",
             def = function(x, y, ...) standardGeneric("plot")
  )
}
