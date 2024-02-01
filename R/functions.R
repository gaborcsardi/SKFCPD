##########################################################################
## SKFCPD construction functions
## 
## SKFCPD Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2023-present Hanmo Li, Yuedong Wang, Mengyang Gu
##							  
##    
##########################################################################

compute_log_lik = function(param, design, response, kernel_type){
  # param: c(log(1/gamma),log(eta))
  
  design = as.matrix(design)
  response = as.matrix(response)
  
  n = nrow(response)
  
  one_matrix = rep(1, n)
  
  fgasp.model_y=fgasp(design, response,kernel_type=kernel_type, have_noise = TRUE)
  fgasp.model_one=fgasp(design, one_matrix,kernel_type=kernel_type, have_noise = TRUE)
  
  log_det_S2_y=Get_log_det_S2_one_dim( param,fgasp.model_y@have_noise,fgasp.model_y@delta_x,
                               fgasp.model_y@output,fgasp.model_y@kernel_type)
  
  log_det_S2_one=Get_log_det_S2_one_dim( param,fgasp.model_one@have_noise,fgasp.model_one@delta_x,
                                 fgasp.model_one@output,fgasp.model_one@kernel_type)
  
  log_det_K = log_det_S2_y[[1]]
  
  y_T_K_inv_y = log_det_S2_y[[2]]
  one_T_K_inv_one = log_det_S2_one[[2]]
  y_T_K_inv_one = sum(log_det_S2_y[[3]] * log_det_S2_one[[3]])
  
  S_2 = (n-1)/2 * log(y_T_K_inv_y - y_T_K_inv_one^2 / one_T_K_inv_one)
  
  log_lik = lgamma((n-1)/2) - 1/2 * log_det_K - 1/2 * log(one_T_K_inv_one) - S_2
  
  return(log_lik)
  
}


get_mu_sigma_hat = function(param, design, response, kernel_type){
  design = as.matrix(design)
  response = as.matrix(response)
  
  n = nrow(response)
  
  one_matrix = rep(1, n)
  
  fgasp.model_y=fgasp(design, response,kernel_type=kernel_type, have_noise = TRUE)
  fgasp.model_one=fgasp(design, one_matrix,kernel_type=kernel_type, have_noise = TRUE)
  
  log_det_S2_y=Get_log_det_S2_one_dim( param,fgasp.model_y@have_noise,fgasp.model_y@delta_x,
                               fgasp.model_y@output,fgasp.model_y@kernel_type)
  
  log_det_S2_one=Get_log_det_S2_one_dim( param,fgasp.model_one@have_noise,fgasp.model_one@delta_x,
                                 fgasp.model_one@output,fgasp.model_one@kernel_type)
  
  y_T_K_inv_y = log_det_S2_y[[2]]
  one_T_K_inv_one = log_det_S2_one[[2]]
  y_T_K_inv_one = sum(log_det_S2_y[[3]] * log_det_S2_one[[3]])
  
  mu_hat = y_T_K_inv_one / one_T_K_inv_one
  sigma_hat = 1/(n-1) * (y_T_K_inv_y - 2 * mu_hat * y_T_K_inv_one + mu_hat^2 * one_T_K_inv_one)
  
  res = list(mu_hat = mu_hat,
             sigma_hat = sigma_hat)
  return(res)
}


Estimate_GP_params <- function(input, output, kernel_type='matern_5_2'){
  
  input = as.matrix(input)
  output = as.matrix(output)
  n_dim = ncol(output)
  
  beta_seq = rep(NA, n_dim)
  eta_seq = rep(NA, n_dim)
  sigma_2_seq = rep(NA, n_dim)
  
  
  for(i in 1:n_dim){
    ini_beta_seq = seq(0.1, 1, 0.2)
    ini_eta_seq = seq(0.02, 0.5, 0.1)
    
    ini_params_set = expand.grid(beta = ini_beta_seq, eta = ini_eta_seq)
    
    final_est_all = NA
    max_log_lik = -Inf
    
    for (param_index in 1:nrow(ini_params_set)) {
      
      ##log inverse range and log noise-variance ratio (nugget) parameters
      param_1 = c(log(ini_params_set$beta[param_index]), log(ini_params_set$eta[param_index]))
      
      est_all <- tryCatch(
        {
          optim(
            param_1,
            compute_log_lik,
            design = input,
            response = output[,i],
            kernel_type = kernel_type,
            method = "L-BFGS-B",
            control = list(fnscale = -1),
            lower = log(c(10^(-2), 10^(-2))),
            upper = log(c(10^(2), 10^(2)))
          )
        },
        error = function(e) {
          NULL  # Return NULL when an error occurs
        }
      )
      
      # Check if the result is not NULL (i.e., no error occurred)
      if (!is.null(est_all)){
        if(est_all$value > max_log_lik){
          max_log_lik = est_all$value
          final_est_all = est_all
        }
      }
    }

    ##estimate variance
    mu_sigma_hat = get_mu_sigma_hat(param = final_est_all$par, 
                              design = input,
                              response = output[,i],
                              kernel_type = kernel_type)
    
    beta_seq[i] = exp(final_est_all$par[1])
    eta_seq[i] = exp(final_est_all$par[2])
    sigma_2_seq[i] = mu_sigma_hat$sigma_hat
    
  }
  
  param_est = new("Estimated_GP_params")
  ##estimated inverse range parameter and nugget
  param_est@beta = beta_seq
  param_est@eta = eta_seq
  param_est@sigma_2 = sigma_2_seq
  return(param_est)
}


SKFCPD <- function(design = NULL, response = NULL, FCPD = NULL, init_params = list(gamma = 1, sigma_2 = 1, eta = 1), train_prop = NULL, kernel_type = "matern_5_2", hazard_vec=100, print_info = TRUE,
                    truncate_at_prev_cp = FALSE){
  
  if( (kernel_type!='matern_5_2')&(kernel_type!='exp') ){
    stop("the current version only support the Matern covariance with smoothness 
         parameter being 2.5 or 0.5 (exponential kernel). \n")
  }
  
  n_obs = nrow(as.matrix(response))
  
  if(is.null(FCPD)){
    if (length(hazard_vec)!=n_obs){
      if(length(hazard_vec)==1){
        hazard_vec = rep(hazard_vec, n_obs)
      } else{
        stop("The hazard_vec needs to be a vector with the length equals to number of observations. \n")
      }
    }
  }
  
  test_start = 1
  
  if (is.null(FCPD)){
    
    if(print_info){
      cat("Starting a new changepoint detection model.\n")
    }
    
    design = as.matrix(design)
    response = as.matrix(response)
    
    n_dim = ncol(response)
    if (n_dim == 1){
      missing_index = which(is.na(response))
    } else{
      miss_loc = function(x){which(is.na(x))}
      identicalValue = function(x,y) if (identical(x,y)) x else FALSE
      missing_index_multi = apply(response, 2, miss_loc, simplify = FALSE)
      missing_index = Reduce(identicalValue,missing_index_multi)
    }
    
    if (is.numeric(missing_index)){
      if(length(missing_index)>0){
        response = response[-missing_index,,drop=FALSE]
        design = design[-missing_index,,drop=FALSE]
        if(print_info){
          cat("Find", length(missing_index), "missing values in the response.\n")
          cat("Skip all the missing values.\n")
        }
      } else{
        if(print_info){
          cat("No missing values found in the response.\n")
        }
      }
    } else{
      stop("Find irregular missing patterns across response dimensions.\n")
    }
    
    n_obs = nrow(design)
    
    if (!is.null(train_prop)){
      if(print_info){
        cat("Starting the training process:\n")
        cat("Use first ", train_prop*100, "% observations for training.\n", sep="")
      }
      
      train_index = seq(1, round(n_obs*train_prop))
      test_index = seq(round(n_obs*train_prop)+1, n_obs)
      
      test_start = test_index[1]
      
      x_train = as.matrix(design[train_index,])
      y_train = as.matrix(response[train_index,])
      
      x_test = as.matrix(design[test_index,])
      y_test = as.matrix(response[test_index,])
      
      hazard_vec = hazard_vec[test_index]
      
      fast_model = Estimate_GP_params(x_train, y_train, kernel_type) 
      
      if(print_info){
        cat("Esimated parameters: \n")
        cat("gamma:", round(1/fast_model@beta,5), "\n")
        cat("eta:", round(fast_model@eta,5), "\n")
        cat("sigma_2:", round(fast_model@sigma_2,5), "\n")
      }
      
      gamma = 1/fast_model@beta
      sigma_2 = fast_model@sigma_2
      eta = fast_model@eta
      
      
      if (length(hazard_vec)!=nrow(x_test)){
        stop("The hazard_vec needs to be a vector with the length equals to number of observations. \n")
      }
      
      # implement FastCPD
      ## detect the change points on rest data
      results <- CPD_DLM(design = x_test,
                         response = y_test,
                         gamma = gamma,
                         model_type = 2,
                         mu = rep(0, n_dim), # doesn't matter when model_type=2
                         sigma_2 = sigma_2,
                         eta = eta,
                         kernel_type = kernel_type,
                         stop_at_first_cp = FALSE,
                         hazard_vec = hazard_vec,
                         truncate_at_prev_cp = truncate_at_prev_cp
      )
      
      if(print_info){
        cat("Fast changepoint detection algorithm completed. \n")
      }
      
      
      run_length = unlist(apply(results$run_length_posterior_mat, 2, which.max))
      cp_location = sort(unique(test_index - run_length+1)) 
      cp_location = cp_location[cp_location!=test_index[1]]
      
      if(print_info){
        cat("Estimated changepoints: ", cp_location, "\n")
      }
      
      
      
      
    } else{
      
      if (is.null(init_params)){
        stop("Please specify the initial parameters: gamma, eta, sigma_2.")
      }
      
      if(print_info){
        cat("Esimated parameters: \n")
        cat("gamma: ", round(init_params$gamma,5), "\n")
        cat("eta: ", round(init_params$eta,5), "\n")
        cat("sigma_2: ", round(init_params$sigma_2,5), "\n")
      }
      
      gamma = init_params$gamma
      sigma_2 = init_params$sigma_2
      eta = init_params$eta
      
      n_dim = ncol(response)
      
      if ((length(gamma)!=n_dim)|(length(sigma_2)!=n_dim)|(length(eta)!=n_dim)){
        stop("The vector length of gamma, eta and sigma_2 don't match with the input dimension.")
      }
      
      
      # implement SKF
      ## detect the change points on rest data
      n_dim = ncol(response)
      results <- CPD_DLM(design = design,
                         response = response,
                         gamma = gamma,
                         model_type = 2,
                         mu = rep(1, n_dim), # doesn't matter when model_type=2
                         sigma_2 = sigma_2,
                         eta = eta,
                         kernel_type = kernel_type,
                         stop_at_first_cp = FALSE,
                         hazard_vec = hazard_vec,
                         truncate_at_prev_cp = truncate_at_prev_cp
      )
      
      run_length = unlist(apply(results$run_length_posterior_mat, 2, which.max))
      cp_location = sort(unique(1:n_obs - run_length+1)) 
      cp_location = cp_location[cp_location!=1]
    }
  } else{
    # if FCPD is not NULL
    
    prev_design = FCPD@design
    pred_response = FCPD@response
    n_dim = ncol(pred_response)
    n_obs = nrow(pred_response)+1
    hazard_vec = c(FCPD@hazard_vec, hazard_vec)
    
    gamma = FCPD@gamma
    sigma_2 = FCPD@sigma_2
    eta = FCPD@eta
    
    if ((length(gamma)!=n_dim)|(length(sigma_2)!=n_dim)|(length(eta)!=n_dim)){
      stop("The vector length of gamma, eta and sigma_2 don't match with the input dimension.")
    }
    
    results = list(KF_params_list = vector("list", n_dim),
                   prev_L_params_list = vector("list", n_dim),
                   # G_W_W0_V_ini = FCPD@G_W_W0_V_ini,
                   run_length_posterior_mat = matrix(NA, n_obs, n_obs),
                   #run_length_joint_mat = matrix(NA, n_obs, n_obs),
                   log_pred_dist_mat = matrix(NA, n_obs, n_obs))
    
    results$log_pred_dist_mat[1:(n_obs-1), 1:(n_obs-1)] = FCPD@log_pred_dist_mat
    results$run_length_posterior_mat[1:(n_obs-1), 1:(n_obs-1)] = FCPD@run_length_posterior_mat
    
    cur_d = as.numeric(abs(design - prev_design[n_obs-1,])) # L1 distance
    
    G_W_W0_V_ini_list = vector("list", n_dim)
    G_W_W0_V_list = vector("list", n_dim)
    
    for(i in 1:n_dim){
      G_W_W0_V_ini_list[[i]] = Construct_G_W_W0_V(d = 1,
                                                  gamma = FCPD@gamma[i],
                                                  eta = FCPD@eta[i],
                                                  kernel_type = kernel_type,
                                                  is_initial = TRUE)
      
      G_W_W0_V_list[[i]] = Construct_G_W_W0_V(d = cur_d,
                                              gamma = FCPD@gamma[i],
                                              eta = FCPD@eta[i],
                                              kernel_type = kernel_type,
                                              is_initial = FALSE)
    }
    
    KF_results_list = vector("list", n_dim)
    
    for (j in 1:n_dim){
      KF_results_list[[j]] = GaSP_CPD_pred_dist_objective_prior_KF_online(KF_params = FCPD@KF_params_list[[j]],
                                                                          # det_params = det_params_list[[j]], 
                                                                          prev_L_params = FCPD@prev_L_params_list[[j]],
                                                                          cur_point = response[j], 
                                                                          d = cur_d,
                                                                          model_type = 2,
                                                                          gamma = gamma[j],
                                                                          mu = 0,
                                                                          sigma_2 = sigma_2[j],
                                                                          eta = eta[j],
                                                                          kernel_type = kernel_type,
                                                                          G_W_W0_V_ini = G_W_W0_V_ini_list[[j]],
                                                                          G_W_W0_V = G_W_W0_V_list[[j]])
      results$KF_params_list[[j]] = KF_results_list[[j]]$KF_params
      # det_params_list[[j]] = KF_results_list[[j]]$det_params
      results$prev_L_params_list[[j]] = KF_results_list[[j]]$cur_L
    }
    
    log_pred_dist_multi_dim = 0
    if(!is.null(KF_results_list[[1]]$log_pred_dist)){
      for(j in 1:n_dim){
        
        log_pred_dist_multi_dim = log_pred_dist_multi_dim + KF_results_list[[j]]$log_pred_dist
      }
      results$log_pred_dist_mat[1:n_obs, n_obs] = log_pred_dist_multi_dim
    }
    
    if (n_obs>=2){
      pred_like = exp(results$log_pred_dist_mat[1:(n_obs-1), n_obs])
      
      # modify p(y_n \mid y_n-1)
      pred_like[length(pred_like)] = 0.2 * pred_like[length(pred_like)]
      
      run_length_joint_vec = rep(0, n_obs)
      # when r_t = r_t-1 + 1
      run_length_joint_vec[1:(n_obs-1)] = (1 - 1/hazard_vec[1:(n_obs-1)]) * pred_like * results$run_length_posterior_mat[(n_obs-1):1, n_obs-1] # results$run_length_joint_mat[1:(n_obs-1), n_obs-1]
      # when r_t = 0
      # run_length_joint_vec[n_obs] = 1/hazard_vec * sum(pred_like * results$run_length_posterior_mat[(n_obs-1):1, n_obs-1])
      run_length_joint_vec[n_obs] =  prod((sigma_2 * (1+eta))^(-1/2)) * 1/hazard_vec[n_obs] * sum(results$run_length_posterior_mat[(n_obs-1):1, n_obs-1])
      # get run length posterior distribution
      results$run_length_posterior_mat[n_obs:1, n_obs] = run_length_joint_vec / sum(run_length_joint_vec)
    }
    
    
    run_length = unlist(apply(results$run_length_posterior_mat, 2, which.max))
    cp_location = sort(unique(1:n_obs - run_length+1)) 
    cp_location = cp_location[cp_location!=1]
    
    design = rbind(FCPD@design, design)
    response = rbind(FCPD@response, response)
    
  }
  
  SKFCPD_results = new("SKFCPD")
  ## data
  SKFCPD_results@design = design
  SKFCPD_results@response = response
  SKFCPD_results@test_start = test_start
  ## parameters
  SKFCPD_results@gamma = gamma
  SKFCPD_results@eta = eta
  SKFCPD_results@sigma_2 = sigma_2
  SKFCPD_results@kernel_type = kernel_type
  SKFCPD_results@hazard_vec = hazard_vec
  ## results
  SKFCPD_results@KF_params_list = results$KF_params_list
  SKFCPD_results@prev_L_params_list = results$prev_L_params_list
  SKFCPD_results@run_length_posterior_mat = results$run_length_posterior_mat
  SKFCPD_results@log_pred_dist_mat = results$log_pred_dist_mat
  SKFCPD_results@cp = cp_location
  
  return(SKFCPD_results)
}
