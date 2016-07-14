#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Sampling from multivariate Gaussian
// [[Rcpp::export]]
arma::rowvec mvrnorm_cpp(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return(arma::trans(mu + Y) * arma::chol(sigma));
}


// Function to update empirical covariance matrix
// [[Rcpp::export]]
arma::mat update_sigma(arma::mat cov_mat, arma::vec residual, int i){
  arma::mat out = (cov_mat * (i-1) + (i+1)/i*residual * residual.t())/i;
  return(out);
}


// Random Walk Metropolis-Hastings MCMC
/*
 * target is a function that returns the log probability density to sample from; it must take a vector as input
 * theta_init is a vector of initial parameter values
 * sigma is the covariance matrix of the Gaussian transition kernel (proposal density)
 * iterations is the number of iterations
 */
// [[Rcpp::export]]
List rw_mcmc(Function target, arma::vec theta_init, arma::mat sigma, int iterations, int info){
  
  arma::vec theta_i = theta_init; //current value of theta
  NumericMatrix theta_samp = NumericMatrix(iterations,theta_init.n_elem); //store the trace of theta
  double acc_rate = 0; //record acceptance rate
  double target_i; //evaluate target at current theta
  target_i = as<double>(target(theta_i));
  
  
  for(int i = 0; i < iterations; i++){
    
    //print diagnostics
    if((i+1) % info == 0){
      Rcout << "Iteration: " << i+1 << ", acceptance rate: " << acc_rate << std::endl;
    }
    
    //sample from the proposal distribution (transition kernel)
    arma::rowvec theta_star;
    theta_star = mvrnorm_cpp(theta_i,sigma);
    
    //evaluate target at proposed theta
    double target_star;
    target_star = as<double>(target(theta_star));
    
    //evaluate MH acceptance kernel
    double A;
    bool acc;
    if(!std::isfinite(target_star)){
      A = -std::numeric_limits<double>::infinity();
      acc = false; //boolean to record acceptance at each i
    } else {
      //compute A (log of the Bayes factor)
      A = target_star - target_i;
      acc = false; //boolean to record acceptance at each i
    }
    
    //calculate MH acceptance probability
    arma::vec armaRand = arma::randu(1);
    double r_num;
    r_num = as<double>(wrap(armaRand(0)));
    if(r_num < exp(A)){
      theta_i = theta_star.t(); //update current value of theta
      target_i = target_star; //update current value of target
      acc = true; //record acceptance
    }
    
    //update acceptance rate
    if(i == 0){
      acc_rate = 1;
    } else {
      acc_rate = acc_rate + (acc - acc_rate) / i;
    }
    
    theta_samp(i,_) = as<NumericVector>(wrap(theta_i)); //record the current value of theta
  }
  
  return(List::create(Named("trace")=theta_samp,Named("acc_rate")=acc_rate));
}


// Random Walk Metropolis-Hastings MCMC with Adaptive Transition Kernel
/*
 * target is a function that returns the log probability density to sample from; it must take a vector as input
 * theta_init is a vector of initial parameter values
 * sigma is the covariance matrix of the Gaussian transition kernel (proposal density)
 * cooling is the cooling factor for scaling the size of sigma
 * adapt_size is the iteration to begin adapting size of sigma
 * adapt_shape is the iteration to begin adapting shape of sigma
 * iterations is the number of iterations
 * info is the frequency of diagnostic information on the chain
 */
// [[Rcpp::export]]
List adapt_mcmc(Function target, arma::vec theta_init, arma::mat sigma, double cooling, int adapt_size, int adapt_shape, int iterations, int info){

  //prepare trace
  arma::vec theta_i = theta_init; //current value of theta
  NumericMatrix theta_samp = NumericMatrix(iterations,theta_init.n_elem); //store the trace of theta
  double acc_rate = 0; //record acceptance rate
  double target_i; //evaluate target at current theta
  target_i = as<double>(target(theta_i));

  //prepare sigma and adaptive routine
  int num_param = theta_init.n_elem;
  int s_dim = sigma.n_rows;
  arma::mat sigma_proposal = sigma; //proposal sigma
  arma::mat sigma_empirical = arma::zeros(s_dim,s_dim); //empirical sigma
  arma::mat sigma_init = sigma_proposal; //initial sigma
  arma::vec theta_mean = theta_i; //empirical mean theta
  arma::cube sigma_trace = arma::zeros(s_dim,s_dim,iterations); //array to return trace of adaptive sigma
  bool adapting_size = false; //boolean to control size adaptation
  bool adapting_shape = false; //boolean to control shape adaptation
  double scale_sd = 1; //scaling factor for size adaptation
  double max_scale_sd = 50; //maximum allowed value for size scaling factor MAYBE TURN INTO PARAMETER LATER
  double scale_multiplier = 1; //scaling multiplier for size adaptation

  //run mcmc loop
  for(int i=0; i < iterations; i++){
    
    //adapt transition kernel size
    if(adapt_size != 0 && adapt_size <= i && (adapt_shape == 0 || acc_rate*i < adapt_shape)){ //adapt size of sigma until enough transitions are accepted
      if(!adapting_size){ //on first iteration of size adaptation print to R console and change the control boolean
        Rcout << "Begin adapting size of sigma at iteration: " << i << std::endl;
        adapting_size = true;
      }
      scale_multiplier = exp(pow(cooling,(i-adapt_size)) * (acc_rate - 0.234));
      scale_sd = scale_sd * scale_multiplier;
      scale_sd = std::min(scale_sd,max_scale_sd);
      arma::mat sigma_new = pow(scale_sd,2)*sigma_init;
      if(!any(sigma_new.diag() < 2E-16)){
        sigma_proposal = sigma_new;
      }
    } else if(adapt_shape != 0 && acc_rate*i >= adapt_shape){ //adapt shape of sigma after enough accepted transitions
      if(!adapting_shape){ //on first iteration of shape adaptation print to R console and change the control boolean
        Rcout << "Begin adapting shape of sigma at iteration: " << i << std::endl;
        adapting_shape = true;
      }
      scale_sd = 2.38/sqrt(num_param);
      sigma_proposal = pow(scale_sd,2)*sigma_empirical;
    }
    
    //print diagnostics
    if((i+1) % info == 0){
      Rcout << "Iteration: " << i+1 << ", acceptance rate: " << acc_rate << std::endl;
    }
    
    //sample from the proposal distribution (transition kernel)
    arma::rowvec theta_star;
    theta_star = mvrnorm_cpp(theta_i,sigma_proposal);
    
    //evaluate target at proposed theta
    double target_star;
    target_star = as<double>(target(theta_star));
    
    //evaluate MH acceptance kernel
    double A;
    bool acc;
    if(!std::isfinite(target_star)){
      A = -std::numeric_limits<double>::infinity();
      acc = false; //boolean to record acceptance at each i
    } else {
      //compute A (log of the Bayes factor)
      A = target_star - target_i;
      acc = false; //boolean to record acceptance at each i
    }
    
    //calculate MH acceptance probability
    arma::vec armaRand = arma::randu(1);
    double r_num;
    r_num = as<double>(wrap(armaRand(0)));
    if(r_num < exp(A)){
      theta_i = theta_star.t(); //update current value of theta
      target_i = target_star; //update current value of target
      acc = true; //record acceptance
    }
    
    //update acceptance rate
    if(i == 0){
      acc_rate = 1;
    } else {
      acc_rate = acc_rate + (acc - acc_rate) / i;
    }
    
    //update empirical covariance matrix (sigma)
    if(i == 0){ //not a great solution ASK HECTOR??
      theta_mean = theta_mean;
      sigma_empirical = sigma_empirical;
    } else {
      arma::vec residual = theta_i - theta_mean;
      theta_mean = theta_mean + (residual/i); //update empirical mean theta
      sigma_empirical = update_sigma(sigma_empirical,residual,i); //update empirical covariance matrix
    }
    
    sigma_trace.slice(i) = sigma_proposal; //record current sigma
    theta_samp(i,_) = as<NumericVector>(wrap(theta_i)); //record the current value of theta
  }

  return(List::create(Named("theta_trace")=theta_samp,Named("sigma_trace")=sigma_trace,Named("acc_rate")=acc_rate));
}



