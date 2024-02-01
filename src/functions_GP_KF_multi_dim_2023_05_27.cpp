
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#define STRICT_R_HEADERS

#include <RcppEigen.h>
#include <Rcpp.h>
#include <cmath>
// #include "ctools.h"



using namespace Rcpp;
using namespace std;
using namespace Eigen;

//#define M_PI 3.14159265359

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//April 21, 2023

// [[Rcpp::export]]
double productPowerMinusHalf(const Rcpp::NumericVector& vec) {
  double product = 1.0;
  
  for (int i = 0; i < vec.size(); ++i) {
    product *= std::pow(vec[i], -0.5);
  }
  
  return product;
}



////Construct_W0_matern_5_2_one_dim 
// [[Rcpp::export]] 
Eigen::MatrixXd Construct_W0_matern_5_2_one_dim(const double lambda){ 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(3,3); 
  
  W0(0,0)=1; 
  W0(0,2)=W0(2,0)=-pow(lambda,2.0)/3.0; 
  W0(1,1)=pow(lambda,2.0)/3.0; 
  W0(2,2)=pow(lambda,4.0); 
  
  return W0; 
} 



////Construct_G_matern_5_2_one_dim 
// [[Rcpp::export]] 
Eigen::MatrixXd Construct_G_matern_5_2_one_dim(double delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  d(0,0)=pow(lambda,2.0)*pow(delta_x,2.0)+2*lambda*delta_x+2; 
  d(1,0)=-pow(lambda,3.0)*pow(delta_x,2.0); 
  d(2,0)=pow(lambda,4.0)*pow(delta_x,2.0)-2*pow(lambda,3.0)*delta_x; 
  d(0,1)=2*(lambda*pow(delta_x,2.0)+delta_x); 
  d(1,1)=-2*(pow(lambda,2.0)*pow(delta_x,2.0)-lambda*delta_x-1); 
  d(2,1)=2*(pow(lambda,3.0)*pow(delta_x,2.0)-3*pow(lambda,2.0)*delta_x); 
  d(0,2)=pow(delta_x,2); 
  d(1,2)=2*delta_x-lambda*pow(delta_x,2.0); 
  d(2,2)=pow(lambda,2.0)*pow(delta_x,2.0)-4*lambda*delta_x+2;     
  d=exp(-lambda*delta_x)/2.0*d;
  
  return d; 
} 

////Construct_W_matern_5_2_one_dim  
// [[Rcpp::export]] 
Eigen::MatrixXd Construct_W_matern_5_2_one_dim(double delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  double  lambda_delta_x;
  double exp_neg_2_lambda_delta_x;
  
  lambda_delta_x=lambda*delta_x;  //close and jump then it is... 
  exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
  
  d(0,0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda,5.0)); 
  d(1, 0)=  d(0, 1)=exp_neg_2_lambda_delta_x*pow(delta_x,4.0)/2.0; 
  d(2, 0)=  d(0, 2)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda,3.0)); 
  d(1, 1)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda,3.0)); 
  d(2, 1)=  d(1, 2)=exp_neg_2_lambda_delta_x*pow(delta_x,2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
  d(2, 2)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda)     ;  
  d=d*(4*pow(lambda,5.0)/3.0); 
  
  return d; 
}

////Construct_G_W_W0_V 
// [[Rcpp::export]] 
List Construct_G_W_W0_V(double d, double gamma, double eta, const String kernel_type,
                        bool is_initial){
  List G_W_W0_V;
  Eigen::MatrixXd G;
  Eigen::MatrixXd W;
  Eigen::MatrixXd W_0;
  double VV;
  double lambda;
  
  if (kernel_type == "exp"){
    
    lambda = 1.0 / gamma;
    G = W = W_0 = Eigen::MatrixXd::Zero(1,1);
    
    if (is_initial){
      W(0,0) = 1.0;
      W_0(0,0) = 1.0;
      VV = eta;
    } else{
      G(0,0) = exp(-1.0 * d / gamma);
      W(0,0) = 1.0 - exp(-2.0 * d / gamma);
      W_0(0,0) = 1.0;
      VV = eta;
    }
  } else if(kernel_type == "matern_5_2"){
    
    lambda = 1.0 * sqrt(5) / gamma;
    G = W = W_0 = Eigen::MatrixXd::Zero(3,3);
    
    if (is_initial){
      W = W_0 = Construct_W0_matern_5_2_one_dim(lambda);
      VV = eta;
    } else{
      G = Construct_G_matern_5_2_one_dim(d, lambda);
      W = Construct_W_matern_5_2_one_dim(d, lambda);
      W_0 = Construct_W0_matern_5_2_one_dim(lambda);
      VV = eta;
    }
  }
  
  G_W_W0_V = List::create(Named("G") = G,
                          Named("W") = W,
                          Named("W_0") = W_0,
                          Named("VV") = VV);
  
  return G_W_W0_V;
  
  
}

////Construct_G_matern_5_2_fastGP 
// [[Rcpp::export]] 
List Construct_G_matern_5_2_fastGP(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(3,3); 
  
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  

  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    int j_GG_1=j_GG+1;    
    d(0,0)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)+2*lambda*delta_x[j_GG]+2; 
    d(1,0)=-pow(lambda,3.0)*pow(delta_x[j_GG],2.0); 
    d(2,0)=pow(lambda,4.0)*pow(delta_x[j_GG],2.0)-2*pow(lambda,3.0)*delta_x[j_GG]; 
    d(0,1)=2*(lambda*pow(delta_x[j_GG],2.0)+delta_x[j_GG]); 
    d(1,1)=-2*(pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-lambda*delta_x[j_GG]-1); 
    d(2,1)=2*(pow(lambda,3.0)*pow(delta_x[j_GG],2.0)-3*pow(lambda,2.0)*delta_x[j_GG]); 
    d(0,2)=pow(delta_x[j_GG],2); 
    d(1,2)=2*delta_x[j_GG]-lambda*pow(delta_x[j_GG],2.0); 
    d(2,2)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-4*lambda*delta_x[j_GG]+2;     
    d=exp(-lambda*delta_x[j_GG])/2.0*d;
    GG[j_GG_1]=d; 
  } 

  return GG; 
} 


////Construct_W_matern_5_2_fastGP  
// [[Rcpp::export]] 
List Construct_W_matern_5_2_fastGP(Eigen::VectorXd delta_x, double lambda, Eigen::MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
 
  double  lambda_delta_x;
  double exp_neg_2_lambda_delta_x;
  int  j_Wi_1;
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    j_Wi_1= j_Wi+1; 
    lambda_delta_x=lambda*delta_x[j_Wi]; 
    exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
    
    d(0,0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda,5.0)); 
    d(1, 0)=  d(0, 1)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],4.0)/2.0; 
    d(2, 0)=  d(0, 2)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda,3.0)); 
    d(1, 1)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda,3.0)); 
    d(2, 1)=  d(1, 2)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
    d(2, 2)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda)     ;  
    d=d*(4*pow(lambda,5.0)/3.0); 
    Wi[j_Wi_1]=d; 
    

  } 
  return Wi; 
}




////Construct_G_exp_fastGP
// [[Rcpp::export]] 
List Construct_G_exp_fastGP(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(1,1); 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    d(0,0)=exp(-delta_x[j_GG]*lambda);
    GG[j_GG+1]=d; 
  }
  
  return GG;
}


////Construct_W_exp_fastGP  
// [[Rcpp::export]] 
List Construct_W_exp_fastGP( Eigen::VectorXd delta_x, double lambda, Eigen::MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
 
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    d(0,0)=1-exp(-2*delta_x[j_Wi]*lambda);
    Wi[j_Wi+1]=d;
  }
  
  return Wi;
}

// [[Rcpp::export]] 
Eigen::MatrixXd Construct_W0_exp_one_dim(const double lambda){ 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(1,1); 
  
  W0(0,0)=1.0; 
  
  return W0; 
} 


////Get_Q_K  
// [[Rcpp::export]] 
List Get_Q_K(const List GG,const List  W,const Eigen::MatrixXd C0,const double VV){ 
  
  int n=GG.size();
  int k=C0.rows();
  
  Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
  Eigen::MatrixXd C=C0;
  
  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd W_matrix;
  
  Eigen::MatrixXd RR;
  
  
  // num_dim list, each is 3(num_obs)\times 3 list 
  for(int t=0;t<n;t++){ 
    GG_matrix=GG[t];
    W_matrix=W[t];
    RR=GG_matrix*C*GG_matrix.transpose()+W_matrix;
    //Q[t]=RR(0,0);
    Q[t]=RR(0,0)+VV;
    K.row(t)=RR.col(0).transpose()/Q[t];
    C=RR-RR.col(0)*RR.row(0)/Q[t];
  }
  
  List return_list;
  return_list.push_back(Q);
  return_list.push_back(K);
  
  return return_list;
}



////Get_log_det_S2_one_dim
// [[Rcpp::export]] 
List Get_log_det_S2_one_dim(const Eigen::VectorXd param,const bool have_noise,const Eigen::VectorXd delta_x,const Eigen::VectorXd output,
                    const  String kernel_type){
  
  int n=output.rows();
  
  double gamma=1.0/exp(param[0]);
  
  double VV=0;
  if(have_noise){
    VV=exp(param[1]);
  }
  
  Eigen::MatrixXd    W0;
  List    GG;
  
  List    W;
  List    Q_K;
  
  
  double lambda=0;
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
    
    W0=Construct_W0_matern_5_2_one_dim(lambda);   
    GG=Construct_G_matern_5_2_fastGP(delta_x,lambda);  
    W=Construct_W_matern_5_2_fastGP(delta_x,lambda,W0);
    
  }else if(kernel_type=="exp"){
    
    lambda=1.0/gamma;
    W0=Construct_W0_exp_one_dim(lambda);  
    GG=Construct_G_exp_fastGP(delta_x,lambda);  
    W=Construct_W_exp_fastGP(delta_x,lambda,W0);
    
  }
  
  Q_K=Get_Q_K(GG,W,W0,VV);
  
  
  Eigen::VectorXd Q=Q_K[0];
  Eigen::MatrixXd K=Q_K[1];
  
  double log_det_R=Q.array().log().sum();
  List return_vec;
  return_vec.push_back(log_det_R);
  
  Eigen::MatrixXd GG_matrix;
  
  Eigen::MatrixXd m=Eigen::VectorXd::Zero(n);
  
  Eigen::VectorXd a;
  
  Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n); 
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    Y_minus_a_1[t]=(output[t]-a[0]);
    
    m=a+K.row(t).transpose()*(output[t]-a[0]);
  }
  
  
  
  double S2=(Y_minus_a_1.array()*Y_minus_a_1.array()/Q.array()).sum(); 
  
  return_vec.push_back(S2);
  
  Eigen::VectorXd L2=(Y_minus_a_1.array()/Q.array().sqrt()); 
  
  return_vec.push_back(L2);
  
  return return_vec;
}

////KF_ini
// [[Rcpp::export]]
List KF_ini(double cur_input, double d, double gamma, double eta, const String kernel_type, List G_W_W0_V) {
  
  List param;
  
  if(kernel_type == "exp"){
    double VV = G_W_W0_V[3]; // eta
    double a_KF = 0.0;
    double R_KF = 1.0;
    
    double f_KF = a_KF;
    double Q_KF = R_KF + VV; // eta = sigma_2_0 / sigma_2
    double Q_KF_inv = 1.0 / Q_KF;    // 1d ##solve(Q_KF[1])
    
    Eigen::MatrixXd m_KF = Eigen::MatrixXd::Zero(1,1);
    Eigen::MatrixXd C_KF = Eigen::MatrixXd::Zero(1,1);
    
    m_KF(0,0) = a_KF + R_KF * Q_KF_inv * (cur_input - f_KF);
    C_KF(0,0) = R_KF - R_KF * Q_KF_inv * R_KF;
    
    param = List::create(Named("f_KF") = f_KF,
                         Named("Q_KF") = Q_KF,
                         Named("m_KF") = m_KF,
                         Named("C_KF") = C_KF);
    
  } else if(kernel_type == "matern_5_2"){
    double VV = G_W_W0_V[3]; // eta
    Eigen::MatrixXd R_KF = G_W_W0_V[2];
    
    double f_KF = 0.0;
    double Q_KF = R_KF.coeff(0,0) + VV; // eta = sigma_2_0 / sigma_2
    double Q_KF_inv = 1.0 / Q_KF;    // 1d ##solve(Q_KF[1])
    Eigen::MatrixXd m_KF = R_KF.col(0) * Q_KF_inv * (cur_input - f_KF); // + a_KF
    Eigen::MatrixXd C_KF = R_KF - R_KF.col(0) * Q_KF_inv * R_KF.row(0);
    
    param = List::create(Named("f_KF") = f_KF,
                         Named("Q_KF") = Q_KF,
                         Named("m_KF") = m_KF,
                         Named("C_KF") = C_KF);
  }
  
  return param;
}

////KF_ini_for_profile_like
// [[Rcpp::export]]
List KF_ini_for_profile_like(double cur_input, double d, double gamma, double eta, const String kernel_type, List G_W_W0_V) {
  
  List KF_ini_L1 = KF_ini(1, d, gamma, eta, kernel_type, G_W_W0_V);
  List KF_ini_LY = KF_ini(cur_input, d, gamma, eta, kernel_type, G_W_W0_V);
  
  List params = List::create(Named("KF_ini_L1") = KF_ini_L1 ,
                             Named("KF_ini_LY") = KF_ini_LY);
  return params;
}


////get_LY_online
// [[Rcpp::export]]
List get_LY_online(double cur_input,List prev_param,double eta, List G_W_W0_V) {
  
  Eigen::MatrixXd G = G_W_W0_V[0];
  Eigen::MatrixXd W = G_W_W0_V[1];
  Eigen::MatrixXd W_0 = G_W_W0_V[2];
  double VV = G_W_W0_V[3];
  
  Eigen::MatrixXd m_KF_prev = prev_param["m_KF"];
  Eigen::MatrixXd C_KF_prev = prev_param["C_KF"];
  
  
  Eigen::MatrixXd a_KF_i = G * m_KF_prev;
  Eigen::MatrixXd R_KF_i = G * C_KF_prev * G.transpose() + W; // *sigma_2
  
  double f_KF_i = a_KF_i.coeff(0,0);
  double Q_KF_i = R_KF_i.coeff(0,0) + VV;
  
  double Q_KF_inv = 1.0 / (Q_KF_i);  //solve
  
  Eigen::MatrixXd m_KF_i = a_KF_i + R_KF_i.col(0) * Q_KF_inv * (cur_input - f_KF_i);
  Eigen::MatrixXd C_KF_i = R_KF_i - R_KF_i.col(0) * Q_KF_inv * R_KF_i.row(0);
  
  List param = List::create(Named("f_KF") = f_KF_i,
                            Named("Q_KF") = Q_KF_i,
                            Named("m_KF") = m_KF_i,
                            Named("C_KF") = C_KF_i);
  
  return param;
}

// [[Rcpp::export]]
List KF_param_update_for_profile_like(double cur_input, double cur_num_obs, List prev_param, double d, double gamma, double eta, const String kernel_type, List G_W_W0_V) {

  
  List KF_ini_L1_prev = prev_param["KF_ini_L1"];
  List KF_ini_LY_prev = prev_param["KF_ini_LY"];
  
  List KF_ini_L1 = get_LY_online(1.0, KF_ini_L1_prev, eta, G_W_W0_V);
  List KF_ini_LY = get_LY_online(cur_input, KF_ini_LY_prev, eta, G_W_W0_V);
  List params = List::create(Named("KF_ini_L1") = KF_ini_L1 ,
                        Named("KF_ini_LY") = KF_ini_LY);
  return params;
  
}


// [[Rcpp::export]]
List get_predictive_dist_KF_objective_prior(double cur_input,
                                            double cur_num_obs,
                                            List   params,
                                            List   prev_L,
                                            double d,
                                            double gamma,
                                            int model_type,
                                            double mu,
                                            double sigma_2,
                                            double eta,
                                            const String kernel_type) {
  
  List KF_ini_L1 = params["KF_ini_L1"];
  List KF_ini_LY = params["KF_ini_LY"];
  
  Eigen::VectorXd L1;
  Eigen::VectorXd LY;
  
  if (prev_L.length() == 0) {
    
    L1 = Eigen::VectorXd(1);
    LY = Eigen::VectorXd(1);
    
    L1 << (1 - as<double>(KF_ini_L1["f_KF"])) / sqrt(as<double>(KF_ini_L1["Q_KF"]));
    LY << (cur_input - as<double>(KF_ini_LY["f_KF"])) / sqrt(as<double>(KF_ini_LY["Q_KF"]));
    
    
  } else{
    
    Eigen::VectorXd L1_prev = prev_L["L1"];
    Eigen::VectorXd LY_prev = prev_L["LY"];
    
    L1 = Eigen::VectorXd(L1_prev.size() + 1);
    LY = Eigen::VectorXd(LY_prev.size() + 1);
    
    L1 << L1_prev, (1 - as<double>(KF_ini_L1["f_KF"])) / sqrt(as<double>(KF_ini_L1["Q_KF"]));
    LY << LY_prev, (cur_input - as<double>(KF_ini_LY["f_KF"])) / sqrt(as<double>(KF_ini_LY["Q_KF"]));
  }
  
  // update L
  List cur_L = List::create(Named("L1") = L1 ,Named("LY") = LY);
  
  int n_L = L1.size();
  List res;
  
  if(n_L==1){
    res = List::create(Named("cur_L") = cur_L);
  }else{
    
    double log_pred_dist;
    
    double y_t_K_inv_y_prev = LY.head(n_L-1).dot(LY.head(n_L-1));
    double y_t_K_inv_y = y_t_K_inv_y_prev + pow(LY(n_L-1), 2); // sum(LY^2)
    
    double H_t_K_inv_H_prev = L1.head(n_L-1).dot(L1.head(n_L-1));
    double H_t_K_inv_H = H_t_K_inv_H_prev + pow(L1(n_L-1), 2); // sum(L1^2)
    
    double y_t_K_inv_H_prev = LY.head(n_L-1).dot(L1.head(n_L-1));
    double y_t_K_inv_H = y_t_K_inv_H_prev + LY(n_L-1) * L1(n_L-1);
    
    double log_K_div = log(as<double>(KF_ini_L1["Q_KF"]));
    if ((model_type == 0) || (model_type == 2)){
      double y_t_Q_y = y_t_K_inv_y - pow(y_t_K_inv_H, 2) / H_t_K_inv_H;
      double y_t_Q_y_prev = y_t_K_inv_y_prev - pow(y_t_K_inv_H_prev, 2) / H_t_K_inv_H_prev;
      
      if (model_type == 0){
        // # unknown mu
        // # known sigma_2
        log_pred_dist = (-1.0 / 2.0) * log(1.0 * sigma_2) - 1.0 / 2.0 * log_K_div - 1.0 / 2.0 * (log(H_t_K_inv_H) - log(H_t_K_inv_H_prev)) - 1.0/(2.0*sigma_2) * (y_t_Q_y - y_t_Q_y_prev);
      } else if (model_type == 2){
        // # unknown mu
        // # unknown sigma_2
        double epsilon = 1e-9;
        if(n_L == 2){
          log_pred_dist = - 1.0 / 2.0 * log_K_div - 1.0/2.0 * log(H_t_K_inv_H / H_t_K_inv_H_prev) - 1.0/2.0 * y_t_Q_y;
        } else if (abs(y_t_Q_y_prev) <= epsilon){
          log_pred_dist = -std::numeric_limits<double>::infinity();
        } else{
          double S_2_minus = - (1.0 * n_L - 1.0) / 2.0 * log(y_t_Q_y) + (1.0 * n_L-2.0)/2.0 * log(y_t_Q_y_prev);
          double log_Gamma_div = std::lgamma((1.0 * n_L-1.0)/2.0) - std::lgamma((1.0 * n_L-2.0)/2.0);
          
          log_pred_dist = 1.0 * log_Gamma_div - 1.0 / 2.0 * log_K_div - 1.0/2.0 * log(H_t_K_inv_H / H_t_K_inv_H_prev) + S_2_minus;
          
        }
      }
    } else if (model_type==1){
      // # known mu
      // # unknown sigma_2
      double S_2_KF = y_t_K_inv_y - 2.0 * mu * y_t_K_inv_H + pow(mu, 2) * H_t_K_inv_H;
      double S_2_KF_prev = y_t_K_inv_y_prev - 2.0 * mu * y_t_K_inv_H_prev + pow(mu, 2) * H_t_K_inv_H_prev;
      
      double S_2_minus = log(S_2_KF/S_2_KF_prev);
      
      double log_Gamma_div = std::lgamma((1.0 * n_L)/2.0) - std::lgamma((1.0 * n_L-1.0)/2.0);
      
      log_pred_dist = log_Gamma_div - 1.0 / 2.0 * log_K_div + S_2_minus;
    }
    
    res = List::create(Named("y_t_K_inv_y") = y_t_K_inv_y ,
                       Named("H_t_K_inv_H") = H_t_K_inv_H,
                       Named("y_t_K_inv_H") = y_t_K_inv_H,
                       Named("log_pred_dist") = log_pred_dist,
                       // Named("det_K") = det_K,
                       Named("cur_L") = cur_L);
  }
  
  return res;
  
  
}


// [[Rcpp::export]]
double get_predictive_dist_direct_objective_prior(Eigen::VectorXd cur_input_seq,
                                                  Eigen::MatrixXd d,
                                                  double gamma,
                                                  double mu,
                                                  double sigma_2,
                                                  double eta) {
  double log_pred_dist;

  int num_obs = cur_input_seq.size();

  Eigen::VectorXd cur_input_seq_prev = cur_input_seq.head(num_obs-1);

  Eigen::VectorXd H = Eigen::MatrixXd::Ones(num_obs, 1);
  Eigen::VectorXd H_prev = Eigen::MatrixXd::Ones(num_obs-1, 1);

  Eigen::MatrixXd Sigma = (-(d/gamma).array().pow(1)).exp().matrix();
  Sigma=Sigma+eta*Eigen::MatrixXd::Identity(num_obs,num_obs);  //not sure
  Eigen::MatrixXd Sigma_prev = Sigma.block(0,0,num_obs-1, num_obs-1);

  LLT<Eigen::MatrixXd> lltOfR(Sigma);             // compute the cholesky decomposition of R called lltofR
  Eigen::MatrixXd Sigma_chol = lltOfR.matrixL();   //retrieve factor L  in the decomposition

  LLT<Eigen::MatrixXd> lltOfRprev(Sigma_prev);             // compute the cholesky decomposition of R called lltofR
  Eigen::MatrixXd Sigma_prev_chol = lltOfRprev.matrixL();   //retrieve factor L  in the decomposition

  Eigen::VectorXd Sigma_inv_H = (Sigma_chol.transpose().triangularView<Upper>().solve(Sigma_chol.triangularView<Lower>().solve(H))).transpose();
  Eigen::VectorXd Sigma_prev_inv_H = (Sigma_prev_chol.transpose().triangularView<Upper>().solve(Sigma_prev_chol.triangularView<Lower>().solve(H_prev))).transpose();

  Eigen::VectorXd Sigma_inv_y = (Sigma_chol.transpose().triangularView<Upper>().solve(Sigma_chol.triangularView<Lower>().solve(cur_input_seq))).transpose();
  Eigen::VectorXd Sigma_prev_inv_y = (Sigma_prev_chol.transpose().triangularView<Upper>().solve(Sigma_prev_chol.triangularView<Lower>().solve(cur_input_seq_prev))).transpose();

  double y_Q_y = cur_input_seq.dot(Sigma_inv_y)  - 1.0 / Sigma_inv_H.dot(H) * pow(Sigma_inv_y.dot(H),2);
  double y_Q_y_prev = cur_input_seq_prev.dot(Sigma_prev_inv_y) - 1.0 / Sigma_prev_inv_H.dot(H_prev) * pow(Sigma_prev_inv_y.dot(H_prev),2);

  double log_det_div = - 1.0/2.0 * (log( Sigma.determinant() ) - log( Sigma_prev.determinant()));
  double log_H_K_inv_H_div = - 1.0/2.0 * (log(Sigma_inv_H.dot(H)) - log(Sigma_prev_inv_H.dot(H_prev)));

  double log_Gamma_div;
  double S_2;

  if (num_obs == 2){
    log_pred_dist = (-1.0 / 2.0) * log(1.0 * sigma_2) + log_det_div - log(Sigma_inv_H.dot(H)) + log(Sigma_prev_inv_H.dot(H_prev)) - 1.0/(2.0*sigma_2) * (y_Q_y - y_Q_y_prev);

  } else{
    log_Gamma_div = std::lgamma((1.0 * num_obs - 1.0)/2.0) - std::lgamma((1.0 * num_obs-2.0)/2.0); // maybe unstable!!!!

    S_2 = - (1.0 * num_obs-1.0)/2.0 * log(y_Q_y) + (1.0 * num_obs-2.0) / 2.0 * log(y_Q_y_prev);

    log_pred_dist = log_Gamma_div + log_det_div + log_H_K_inv_H_div + S_2;
  }

  return log_pred_dist;

}

// [[Rcpp::export]]
List GaSP_CPD_pred_dist_objective_prior_KF_online(List KF_params,
                                                  List prev_L_params,
                                                  double cur_point,
                                                  double d,
                                                  double gamma,
                                                  int model_type,
                                                  double mu,
                                                  double sigma_2,
                                                  double eta,
                                                  const String kernel_type,
                                                  List G_W_W0_V_ini,
                                                  List G_W_W0_V){
  
  
  // update the params for X_n
  List updated_KF_params = KF_ini_for_profile_like(cur_point,d,gamma,eta,kernel_type, G_W_W0_V_ini);
  
  // compute profile likelihood by KF parameters
  List pre_res;
  pre_res = get_predictive_dist_KF_objective_prior(cur_point,1,updated_KF_params,List(), d,gamma, model_type,mu,sigma_2,eta,kernel_type);
  
  
  List cur_L_prev = as<List>(pre_res["cur_L"]);
  
  if (KF_params.length() == 0) {
    
    List updated_KF_param_seq(1);
    List cur_L_prev_seq(1);
    List results;
    
    updated_KF_param_seq[0] = updated_KF_params;
    cur_L_prev_seq[0] = cur_L_prev;
    
    results = List::create(Named("KF_params") = updated_KF_param_seq ,
                           Named("cur_L") = cur_L_prev_seq);
                          
    return results;
  }
  
  int num_prev_point = KF_params.length();
  
  List res_KF_params(num_prev_point + 1);
  List res_L_params(num_prev_point + 1);
  NumericVector res_log_pred_dist(num_prev_point + 1);

  res_KF_params[num_prev_point] = updated_KF_params;
  res_L_params[num_prev_point] = cur_L_prev;
  
  for (int index=0; index <= (num_prev_point-1); ++index) {
    // update parameters for 1:n, 2:n, ..., (n-1):n
    List cur_KF_param = as<List>(KF_params[index]);
    List cur_L = as<List>(prev_L_params[index]);

    
    List updated_KF_params = KF_param_update_for_profile_like(cur_point,num_prev_point - index + 1,cur_KF_param,d,gamma,eta, kernel_type, G_W_W0_V);
    // compute predictive dist by KF parameters
    pre_res = get_predictive_dist_KF_objective_prior(cur_point,num_prev_point - index + 1,updated_KF_params,cur_L,d,gamma,model_type,mu, sigma_2,eta, kernel_type);
    
    res_L_params[index] = as<List>(pre_res["cur_L"]);
    res_KF_params[index] = updated_KF_params;
    res_log_pred_dist[index] = as<double>(pre_res["log_pred_dist"]);
    
  }
  
  
  List results = List::create(Named("KF_params") = res_KF_params ,
                              Named("cur_L") = res_L_params,
                              // Named("det_params") = res_det_params,
                              Named("log_pred_dist") = res_log_pred_dist);
  
  return results;
  
}

// [[Rcpp::export]]
Eigen::VectorXd GaSP_CPD_pred_dist_objective_prior_direct_online(
    Eigen::VectorXd cur_seq,
    List d,
    double gamma,
    double eta,
    double mu,
    double sigma_2) {

  int n = cur_seq.size();

  Eigen::VectorXd res_log_pred_dist(n);

  for (int i=1; i <= (n-1); ++i) {
    // update parameters for 1:n, 2:n, ..., (n-1):n

    res_log_pred_dist(i-1) = get_predictive_dist_direct_objective_prior(cur_seq.tail(n-i+1),d[i-1],gamma,mu,sigma_2,eta);
  }

  return res_log_pred_dist;
}


// [[Rcpp::export]]
List CPD_DLM(NumericVector design,
            Eigen::MatrixXd response,
            NumericVector gamma,
            int model_type,
            NumericVector mu,
            NumericVector sigma_2,
            NumericVector eta,
            const String kernel_type,
            bool stop_at_first_cp,
            Eigen::VectorXd hazard_vec,
            bool truncate_at_prev_cp){

  int num_obs = response.rows();
  int num_dim = response.cols();
  int num_design = design.length();
  
  if (num_obs != num_design){
    Rcpp::stop("ERROR: number of observations mismatch between design and response.");
  }
  
  List KF_params_list(num_dim);
  List prev_L_params_list(num_dim);
  List output;
  List G_W_W0_V_ini_list(num_dim);
  List G_W_W0_V_list(num_dim);
  
  Eigen::MatrixXd::Index maxIndex;
  
  Eigen::VectorXd cp_location = Eigen::VectorXd::Zero(num_obs);
  cp_location = (cp_location.array() + 1.0).matrix();
  
  Eigen::VectorXd max_runlength_location = Eigen::VectorXd::Zero(num_obs);
  max_runlength_location = (max_runlength_location.array()).matrix();
  
  Eigen::MatrixXd pred_like_all;
  Eigen::MatrixXd pred_like;
  Eigen::MatrixXd joint_dist;
  Eigen::MatrixXd product;
  
  Eigen::MatrixXd log_pred_dist_mat = Eigen::MatrixXd::Zero(num_obs, num_obs);
  Eigen::MatrixXd run_length_posterior_mat = Eigen::MatrixXd::Zero(num_obs, num_obs);
  run_length_posterior_mat(0,0) = 1.0;
  run_length_posterior_mat(1,1) = 1.0;
  run_length_posterior_mat(2,2) = 1.0;
  
  
  for(int i = 0; i<num_dim; ++i){
    // double d, double gamma, double eta, const String kernel_type,
    // bool is_initial
    G_W_W0_V_ini_list[i] = Construct_G_W_W0_V(1.0, // double d, doesn't matter in initial state
                                              gamma[i],
                                                   eta[i],
                                                      kernel_type,
                                                      TRUE);

  }
  
  for(int i = 0; i<num_obs; ++i){
    
    double d = 1.0; // doesn't matter in initial state
    if (i > 0){
      d = abs(design[i] - design[i-1]);
    }

    List KF_results_list(num_dim);
    
    for(int j = 0; j<num_dim; ++j){
      
      G_W_W0_V_list[j] = Construct_G_W_W0_V(d,
                                            gamma[j],
                                                 eta[j],
                                                    kernel_type,
                                                    FALSE);
      
      KF_results_list[j] = GaSP_CPD_pred_dist_objective_prior_KF_online(as<List>(KF_params_list[j]), // List KF_params,
                                                                       as<List>(prev_L_params_list[j]), // List prev_L_params,
                                                                       response(i, j),     // double cur_point,
                                                                        d,    // double d,
                                                                        gamma[j],    // double gamma,
                                                                             2,    // int model_type,
                                                                             mu[j],    // double mu,
                                                                               sigma_2[j],    // double sigma_2,
                                                                                      eta[j],    // double eta,
                                                                                         kernel_type,    // const String kernel_type,
                                                                                         as<List>(G_W_W0_V_ini_list[j]),    // List G_W_W0_V_ini,
                                                                                         as<List>(G_W_W0_V_list[j]));   // List G_W_W0_V); 
      KF_params_list[j] = as<List>(KF_results_list[j])["KF_params"];
      prev_L_params_list[j] = as<List>(KF_results_list[j])["cur_L"];
    }
    
    
    
    if(as<List>(KF_results_list[0]).containsElementNamed("log_pred_dist")){
      Eigen::VectorXd log_pred_dist_multi_dim = as<Eigen::Map<Eigen::VectorXd>>(as<List>(KF_results_list[0])["log_pred_dist"]);
      if(num_dim > 1){
        for(int j = 1; j<num_dim; ++j){
          log_pred_dist_multi_dim = log_pred_dist_multi_dim + as<Eigen::Map<Eigen::VectorXd>>(as<List>(KF_results_list[j])["log_pred_dist"]);
        }
      }

      
      
      log_pred_dist_mat.col(i).head(i+1) << log_pred_dist_multi_dim;
    }
    
    if (i>=1){
      
      int prev_max_index = i;
      
      if (truncate_at_prev_cp){

        prev_max_index = max_runlength_location(i-1)+1;


        // [p(y_n | y_{\hat{C}_{n-1}:n-1}), ..., p(y_n | y_{n-1})]
        pred_like = log_pred_dist_mat.col(i).middleRows(i-prev_max_index, prev_max_index).array().exp().matrix();
        
        // modify p(y_n | y_n-1)
        pred_like.col(0).tail(1) = 0.2 * pred_like.col(0).tail(1);
        
        // [p(y_{1:n-1}, C_{n-1}=\hat{C}_{n-1}), ..., p(y_{1:n-1}, C_{n-1}=n-3)]
        joint_dist = run_length_posterior_mat.col(i-1).middleRows(0, prev_max_index).reverse();
        
        product =  pred_like.cwiseProduct(joint_dist);
        
        // example:
        // when i=6, we have n=i+1=7 observations
        // max_runlength_location(5)=2, meaning r_6=3 and C_6=4. Thus t=4 is a detected changepoint
        // we have
        // prev_max_index = 3
        // pred_like_all = [p(y_7 | y_{4:6}), p(y_7 | y_{5:6}), p(y_7 | y_{6})]
        // pred_like = [p(y_7 | y_{4:6}), p(y_7 | y_{5:6}), p(y_7 | y_{6})]
        // joint_dist = [p(y_6, C_6=4), p(y_6, C_6=5), p(y_6, C_6=6)] = rev( [p(y_6, r_6=1), p(y_6, r_6=2), p(y_6, r_6=3)] )
      } else{
        
        // [p(y_n | y_{1:n-1}), ..., p(y_n | y_{n-1})]
        pred_like = log_pred_dist_mat.col(i).head(i).array().exp().matrix();
        
        // modify p(y_n | y_n-1)
        pred_like.col(0).tail(1) = 0.2 * pred_like.col(0).tail(1);
        
        // [p(y_{1:n-1}, C_{n-1}=1), ..., p(y_{1:n-1}, C_{n-1}=n-1)]
        joint_dist = run_length_posterior_mat.col(i-1).middleRows(0, i).reverse();
          
        // product =  p(y_n | y_{k:n-1}) * p(y_{1:n-1}, C_{n-1}=k) for k=1, ..., n-3
        product =  pred_like.cwiseProduct(joint_dist);
        
        // example:
        // when i=6, we have n=i+1=7 observations
        // pred_like_all = [p(y_7 | y_{1:6}), ..., p(y_7 | y_{5:6}), p(y_7 | y_{6})]
        // pred_like = [p(y_7 | y_{1:6}), ..., p(y_7 | y_{6})]
        // joint_dist = [p(y_1:6, C_6=1), ..., p(y_1:6, C_6=4)] = rev( [p(y_1:6, r_6=1), ..., p(y_1:6, r_6=6)] )
      }


      // [p(y_{1:n}, C_n=1), ..., p(y_{1:n}, C_n=n)] = [p(y_{1:n}, r_n=n), ..., p(y_{1:n}, r_n=1)]
      Eigen::VectorXd run_length_joint_vec = Eigen::VectorXd::Zero(prev_max_index+1);
      run_length_joint_vec = (run_length_joint_vec.array()).matrix();
      Eigen::VectorXd one_vec = Eigen::VectorXd::Ones(prev_max_index);
      
      // when r_t = r_t-1 + 1
      run_length_joint_vec.head(prev_max_index) = (one_vec - one_vec.cwiseQuotient(hazard_vec.middleRows(i-prev_max_index, prev_max_index))).cwiseProduct(product); // r_t from n to 3
      // run_length_joint_vec.head(prev_max_index) = (one_vec - one_vec.cwiseQuotient(hazard_vec.segment(i-prev_max_index, prev_max_index))).cwiseProduct(product);
      // when r_t = 0
      run_length_joint_vec(prev_max_index) = productPowerMinusHalf(sigma_2 * (1.0 + eta))  * 1.0/hazard_vec(i) * joint_dist.sum(); // r_t = 2
      
      // get run length posterior distribution
      run_length_posterior_mat.col(i).middleRows(0, prev_max_index+1) << (run_length_joint_vec / run_length_joint_vec.sum()).reverse();
      //middleRows(i, q): Block containing the q rows starting from i
      
      run_length_posterior_mat.col(i).middleRows(0, prev_max_index+1).maxCoeff( &maxIndex );
      
      max_runlength_location(i) = maxIndex;
      cp_location(i) = i + 1 - (maxIndex);
      
      // Example:
      // when i=6, we have n=i+1=7 observations
      // prev_max_index = i = 6
      // run_length_posterior_mat = [p(y_{1:7}, r_7=3), p(y_{1:7}, r_7=4), p(y_{1:7}, r_7=5), p(y_{1:7}, r_7=6), p(y_{1:7}, r_7=7)]
      // suppose p(y_{1:7}, r_7=3) is the largest value in run_length_posterior_mat
      // we have maxIndex = 0
      // max_runlength_location = maxIndex + 2 = 2, meaning the index of maximum value in max_runlength_matrix is 2
      // cp_location = 6 + 1 - (maxIndex + 2) = 5, meaning t=5 is a changepoint, which aligns with the fact that p(y_{1:7}, r_7=3) is the largest value
    }
    
  }
  
  output = List::create(Named("run_length_posterior_mat") = run_length_posterior_mat,
                        //Named("run_length_joint_mat") = run_length_joint_mat,
                        Named("log_pred_dist_mat") = log_pred_dist_mat,
                        Named("cp") = cp_location,
                        Named("max_runlength_index") = max_runlength_location,
                        Named("pred_like") = pred_like,
                        Named("KF_params_list") = KF_params_list,
                        Named("prev_L_params_list") = prev_L_params_list,
                        Named("G_W_W0_V_ini_list") = G_W_W0_V_ini_list);
  
  return output;
  
  
}

