#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


/*
* The function adaptMCMC is an adaptive Metropolis algorithm (probability.ca/jeff/ftpdir/adaptex.pdf)
* It does not do adaptive rejection yet, but usually adapting the covariance matrix of proposal distribution
* tends to give better results than simply adjusting rejection.
* See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance for how the covariance matrix is estimated
* 
* You will need to install Rcpp and RcppArmadillo packages and have a working C++ compiler to use it
* 
* The function takes the following arguments:
* target: this is a R function that must return the log likelihood of the target distribution we want to sample from (ie; log(likelihood) + log(prior))
* init_theta: a vector of initial parameter values (starting point of chain)
* covmat: covariance matrix of the proposal distribution (it is ok to just use identity matrix, because it will be adapted from data anyway)
* n_iterations: integer value, how long to run the chain
* adapt_size_start: number of accepted jumps after which to begin adapting size of proposal covariance matrix
* adapt_shape_start: number of accepted jumps after which to begin adpating shape of proposal covariance matrix (should set higher than adapt_size_start by 1.5X or 2X at least)
* info: print information on chain every X iterations
* adapt_size_cooling: cooling value for scaling size of covariance matrix (default is 0.99, must set between 0 and 1; usually dont need to change)
* 
* note: if you set adapt_size_start and adapt_shape_start to 0, it will run normal random walk MCMC without adaptive routine
*/

// MCMC utility functions
// [[Rcpp::export]]
arma::vec mvrnorm_samp(arma::vec mu, arma::mat sigma) {
  arma::rowvec Y = rnorm(sigma.n_cols,0,1);
  arma::rowvec out = mu.t() + Y * arma::chol(sigma);
  return(out.t());
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols; //how many columns in sigma (ie; how many dimension)
  arma::mat Y = arma::randn(n, ncols); // n draws from N(0,1) 
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat update_sigma(arma::mat sigma, arma::vec residual, double i){
  arma::mat out = (sigma * (i-1) + (i-1) / i * residual * trans(residual)) / i;
  return(out);
}

// [[Rcpp::export]]
double mvrnorm_pdf(arma::vec x, arma::vec mu, arma::mat sigma){
  
  //define constants
  int k = x.size();
  double twoPi = 2 * M_PI;
  double out;
  
  double rootTerm;
  rootTerm = 1 / sqrt(pow(twoPi,k) * det(sigma));
  arma::mat AexpTerm;
  AexpTerm = exp(-0.5 * arma::trans(x - mu) * inv(sigma) * (x - mu));
  double expTerm;
  expTerm = AexpTerm(0,0);
  
  out = rootTerm * expTerm;
  return(log(out));
}


/*
* Random Walk MCMC with Adaptive Proposal Kernel
* adaptMCMC returns an R list:
* theta_trace a matrix recording the trace of the sampled parameters
* sigma_empirical is an array recording the trace of the estimated covariance matrix of the target
* acceptance_rate is a numeric that is the acceptance rate at the end of the MCMC run
*/
// [[Rcpp::export]]
List adaptMCMC(Function target, arma::vec init_theta, arma::mat covmat, int iterations, int adapt_size_start, double acceptance_rate_weight, int acceptance_window, int adapt_shape_start, int info, double adapt_size_cooling = 0.99, double max_scaling_sd = 50.0){
  
  arma::vec theta_current = init_theta;
  arma:: vec theta_propose = init_theta;
  arma::mat covmat_proposal = covmat;
  
  arma::mat covmat_proposal_init = covmat_proposal;
  bool adapting_size = false;
  bool adapting_shape = false;
  
  arma::mat theta_trace = arma::zeros(iterations,init_theta.n_elem);
  arma::cube sigma_trace = arma::zeros(covmat.n_rows,covmat.n_cols,iterations);
  
  //evaluate target distribution at theta_init
  double target_theta_current; //evaluate target at current theta
  target_theta_current = as<double>(wrap(target(theta_current)));
  
  double acceptance_rate = 0.0;
  arma::vec acceptance_series;
  
  double scaling_sd = 1.0;
  double scaling_multiplier = 1.0;
  arma::mat covmat_empirical = arma::zeros(covmat.n_rows,covmat.n_cols);
  arma::vec theta_mean = theta_current;
  
  //main mcmc loop
  for(int i=1; i<=iterations; i++){
    
    //adaptive routine
    if(adapt_size_start != 0 && i >= adapt_size_start && (adapt_shape_start == 0 || acceptance_rate*i < adapt_shape_start)){
      if(!adapting_size){
        Rcout << "Begin adapting size of sigma at iter: " << i << std::endl;
        adapting_size = true;
      }
      //adapt size of sigma until we get enough accepted jumps
      scaling_multiplier = exp(pow(adapt_size_cooling,i-adapt_size_start) * (acceptance_rate - 0.234));
      scaling_sd = scaling_sd * scaling_multiplier;
      scaling_sd = std::min(scaling_sd,max_scaling_sd);
      //only scale if sigma not reduced to zero
      arma::mat covmat_proposal_new = pow(scaling_sd,2) * covmat_proposal_init;
      if(!any(covmat_proposal_new.diag() < 2E-16)){
        covmat_proposal = covmat_proposal_new;
      }
    } else if(adapt_shape_start != 0 && acceptance_rate*i >= adapt_shape_start){
      if(!adapting_shape){
        Rcout << "Begin adapting shape of sigma at iter: " << i << std::endl;
        adapting_shape = true;
      }
      //adapt shape of sigma using optimal scaling  factor for multivariate target distributions
      scaling_sd = 2.38/sqrt(init_theta.n_elem);
      covmat_proposal = pow(scaling_sd,2) * covmat_empirical;
    }
    
    //print chain diagnostics
    if(i % info == 0){
      // Rcout << "At iter: " << i << ", acceptance rate is: " << acceptance_rate << std::endl;
      Rcout << "At iter: " << i << ", acceptance rate is: " << acceptance_rate << ", scaling_sd: " << scaling_sd << ", scaling_multiplier: " << scaling_multiplier << std::endl;
      
    }
    
    //propose new theta
    theta_propose = mvrnorm_samp(theta_current,covmat_proposal);
    
    //evaluate target distribution at proposed theta
    double target_theta_propose;
    target_theta_propose = as<double>(wrap(target(theta_propose)));
    bool is_accepted;
    double log_acceptance;
    
    //if posterior is 0 immediately reject
    if(!std::isfinite(target_theta_propose)){
      is_accepted = false;
      log_acceptance = -std::numeric_limits<double>::infinity();
    } else {
      //compute Metropolis-Hastings ratio (acceptance probability)
      log_acceptance = target_theta_propose - target_theta_current;
      log_acceptance = log_acceptance + mvrnorm_pdf(theta_current,theta_propose,covmat_proposal);
      log_acceptance = log_acceptance - mvrnorm_pdf(theta_propose,theta_current,covmat_proposal);
    }
    
    //evaluate acceptance probability
    double A = R::runif(0,1);
    if(log(A) < log_acceptance){
      //accept proposed parameter set
      is_accepted = true;
      theta_current = theta_propose;
      target_theta_current = target_theta_propose;
    } else {
      is_accepted = false;
    }
    
    //store trace of MCMC
    theta_trace.row(i-1) = theta_current.t();
    
    //update acceptance rate
    if(i == 1){
      acceptance_rate = is_accepted;
    } else {
      if(acceptance_rate_weight == 0){
        if(acceptance_window == 0){
          acceptance_rate = acceptance_rate + (is_accepted - acceptance_rate) / i;
        } else {
          arma::vec is_accepted_vec(1);
          is_accepted_vec(0) = is_accepted;
          acceptance_series = arma::join_cols<arma::mat>(is_accepted_vec,acceptance_series);
          if(acceptance_series.n_elem > acceptance_window){
            int series_length = acceptance_series.n_elem;
            acceptance_series.resize(series_length-1);
          }
          acceptance_rate = arma::mean(acceptance_series);
        }
      } else {
        acceptance_rate = acceptance_rate * (1 - acceptance_rate_weight) + is_accepted * acceptance_rate_weight;
      }
    }
    
    //update empirical covariance matrix (estimation of sigma)
    arma::vec residual = theta_current - theta_mean;
    covmat_empirical = update_sigma(covmat_empirical,residual,i);
    theta_mean = theta_mean + residual/i;
    
    sigma_trace.slice(i-1) = covmat_empirical;
    
  }
  
  return(List::create(Named("theta_trace")=theta_trace,Named("sigma_trace")=sigma_trace,Named("acceptance_rate")=acceptance_rate));
}


/*
* Simple Random Walk MCMC with Adaptive Proposal Kernel 
* adaptMCMC returns an R list:
* theta_trace a matrix recording the trace of the sampled parameters
* sigma_empirical is an array recording the trace of the estimated covariance matrix of the target
* acceptance_rate is a numeric that is the acceptance rate at the end of the MCMC run
*/
//[[Rcpp::export]]
List adaptMCMC_simple(Function target, arma::vec init_theta, arma::mat covmat, int n_iterations, int adapt_size_start, int adapt_shape_start, int info, bool verbose, double adapt_size_cooling = 0.99){
  
  arma::vec theta_current = init_theta;
  arma::vec theta_propose = init_theta;
  arma::mat covmat_proposal = covmat;
  
  arma::mat covmat_proposal_init = covmat_proposal;
  bool adapting_size = false;
  bool adapting_shape = false;
  
  double target_theta_current = as<double>(wrap(target(theta_current)));
  
  arma::mat theta_trace = arma::zeros(n_iterations,init_theta.n_elem);
  arma::cube sigma_empirical = arma::zeros(covmat.n_rows,covmat.n_cols,n_iterations);
  
  //prepare extra data for verbose output
  arma::vec acceptance_trace;
  arma::vec target_theta_current_trace;
  arma::vec scaling_sd_trace;
  arma::cube covmat_proposal_trace;
  arma::mat residual_trace;
  arma::mat theta_mean_trace;
  
  if(verbose){
    acceptance_trace = arma::zeros(n_iterations);
    target_theta_current_trace = arma::zeros(n_iterations);
    scaling_sd_trace = arma::zeros(n_iterations);
    covmat_proposal_trace = arma::zeros(covmat.n_rows,covmat.n_cols,n_iterations);
    residual_trace = arma::zeros(n_iterations,init_theta.n_elem);
    theta_mean_trace = arma::zeros(n_iterations,init_theta.n_elem);
  }
  
  double acceptance_rate = 0.0;
  double scaling_sd = 1.0;
  
  arma::mat covmat_empirical = arma::zeros(covmat.n_rows,covmat.n_cols);
  arma::vec theta_mean = theta_current;
  
  for(int i=1; i<=n_iterations; i++){
    
    //adaptive mcmc routine
    if(adapt_size_start != 0 && i >= adapt_size_start && acceptance_rate * i < adapt_shape_start){
      if(!adapting_size){
        Rcout << "Begin adapting size of sigma at iter: " << i << std::endl;
        adapting_size = true;
      }
      scaling_sd = scaling_sd * exp(pow(adapt_size_cooling,i - adapt_size_start) * (acceptance_rate - 0.234));
      covmat_proposal = pow(scaling_sd,2) * covmat_proposal_init;
    } else if(adapt_shape_start != 0 && acceptance_rate * i >= adapt_shape_start){
      if(!adapting_shape){
        Rcout << "Begin adapting shape of sigma at iter: " << i << std::endl;
        adapting_shape = true;
      }
      covmat_proposal = pow(2.38,2)/init_theta.n_elem * covmat_empirical;
    }
    
    //print chain info
    if(i % info == 0){
      // Rcout << "At iter: " << i << ", acceptance rate is: " << acceptance_rate << std::endl;
      Rcout << "At iter: " << i << ", acceptance rate is: " << acceptance_rate << ", scaling_sd: " << scaling_sd << std::endl;
    }
    
    theta_propose = mvrnorm_samp(theta_current,covmat_proposal);
    double target_theta_propose;
    target_theta_propose = as<double>(wrap(target(theta_propose)));
    
    bool is_accepted;
    double log_acceptance;
    
    if(!std::isfinite(target_theta_propose)){
      log_acceptance = -std::numeric_limits<double>::infinity();
      is_accepted = false;
    } else {
      log_acceptance = target_theta_propose - target_theta_current;
      log_acceptance = log_acceptance + mvrnorm_pdf(theta_current,theta_propose,covmat_proposal);
      log_acceptance = log_acceptance - mvrnorm_pdf(theta_propose,theta_current,covmat_proposal);
    }
    
    double A = R::runif(0,1);
    if(log(A) < log_acceptance){
      theta_current = theta_propose;
      target_theta_current = target_theta_propose;
      is_accepted = true;
    } else {
      is_accepted = false;
    }
    
    theta_trace.row(i-1) = theta_current.t();
    
    acceptance_rate = acceptance_rate + (is_accepted - acceptance_rate)/i;
    
    arma::vec residual = theta_current - theta_mean;
    covmat_empirical = update_sigma(covmat_empirical,residual,i);
    theta_mean = theta_mean + residual/i;
    
    sigma_empirical.slice(i-1) = covmat_empirical;
    
    //extra output for debugging
    if(verbose){
      acceptance_trace(i-1) = acceptance_rate;
      scaling_sd_trace(i-1) = scaling_sd;
      target_theta_current_trace(i-1) = target_theta_current;
      covmat_proposal_trace.slice(i-1) = covmat_proposal;
      residual_trace.row(i-1) = residual.t();
      theta_mean_trace.row(i-1) = theta_mean.t();
    }

  }
  
  if(verbose){
    return(List::create(Named("theta_trace")=theta_trace,Named("sigma_empirical")=sigma_empirical,Named("acceptance_rate")=acceptance_rate,Named("acceptance_trace")=acceptance_trace,Named("target_theta_current_trace")=target_theta_current_trace,Named("covmat_proposal_trace")=covmat_proposal_trace,Named("residual_trace")=residual_trace,Named("theta_mean_trace")=theta_mean_trace,Named("scaling_sd_trace")=scaling_sd_trace));
  } else {
    return(List::create(Named("theta_trace")=theta_trace,Named("sigma_empirical")=sigma_empirical,Named("acceptance_rate")=acceptance_rate));
  }
}