#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// // Sampling from multivariate Gaussian
// // [[Rcpp::export]]
// arma::rowvec mvrnorm_cpp(arma::vec mu, arma::mat sigma) {
//   arma::vec Y = arma::randn(sigma.n_cols);
//   return(arma::trans(mu + Y) * arma::chol(sigma));
// }


// Sampling from multivariate Gaussian
// [[Rcpp::export]]
arma::vec mvrnorm_cpp(arma::vec mu, arma::mat sigma) {
  arma::vec Y = arma::randn(sigma.n_cols);
  arma::rowvec out = arma::trans(mu + Y) * arma::chol(sigma);
  return(out.t());
}


// Function to update empirical covariance matrix
// [[Rcpp::export]]
arma::mat update_sigma(arma::mat sigma, arma::vec residual, int i){
  arma::mat out = (sigma * (i-1) + (i+1)/i*residual * residual.t())/i;
  return(out);
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
List adapt_mcmc(Function target, arma::vec theta_init, arma::mat sigma, double cooling, int adapt_size, int adapt_shape, int iterations, int info, double max_scaling_sd = 50){
  
  //initialize theta
  arma::vec theta_current = theta_init;
  arma::vec theta_propose = theta_init;
  
  //initialize sigma
  arma::mat sigma_proposal = sigma;
  arma::mat sigma_proposal_init = sigma; 
  
  //adaptive sigma routine
  bool adapting_size = false;
  bool adapting_shape = false;
  
  //evaluate target distribution at theta_init
  double target_theta_current; //evaluate target at current theta
  target_theta_current = as<double>(wrap(target(theta_current)));
  
  //store MCMC trace information
  arma::mat theta_samp = arma::zeros(iterations,theta_init.n_elem); //store the trace of theta
  arma::cube sigma_trace = arma::zeros(theta_init.n_elem,theta_init.n_elem,iterations); //array to return trace of adaptive sigma

  double acceptance_rate = 0.0; //acceptance rate
  double scaling_sd = 1.0; //scaling factor for covariance matrix size
  double scaling_multiplier = 1.0; //scaling multiplier
  
  arma::mat sigma_empirical = arma::zeros(sigma.n_rows,sigma.n_cols); //empirical covariance matrix
  arma::vec theta_mean = theta_current; //empirical mean vector
  
  //main MCMC loop
  for(int i=0; i<iterations; i++){
    
    //adapt sigma
    if(i >= adapt_size && acceptance_rate*i < adapt_shape){ //begin adapting size; continue until shape adaptation begins
      if(!adapting_size){ //print message when size adaptation begins
        Rcout << "Begin adapting size of sigma at iteration: " << i << std::endl;
        adapting_size = true;
      }
      //adapt size of sigma until we get enough accepted jumps
      scaling_multiplier = exp(pow(cooling,i*adapt_size) * (acceptance_rate - 0.234));
      scaling_sd = scaling_sd * scaling_multiplier;
      scaling_sd = std::min(scaling_sd,max_scaling_sd);
      //only scale if it doesn't reduce sigma to zero
      arma::mat sigma_proposal_new = pow(scaling_sd,2)*sigma_proposal_init;
      if(!any(sigma_proposal_new.diag() < 2E-16)){
        sigma_proposal = sigma_proposal_new;
      }
    } else if(acceptance_rate*i >= adapt_shape){
      if(!adapting_shape){ //print message when shape adaptation begins
        Rcout << "Begin adapting shape of sigma at iteration: " << i << std::endl;
        adapting_shape = true;
      }
      //adapt shape of sigma using optimal scaling factor for multivariate target distribution
      scaling_sd = 2.38/sqrt(theta_init.n_elem);
      sigma_proposal = pow(scaling_sd,2) * sigma_empirical;
    }
    
    //print diagnostics
    if((i+1) % info == 0){
      Rcout << "Iteration: " << i+1 << ", acceptance rate: " << acceptance_rate << std::endl;
    }
    
    //propose another parameter set
    theta_propose = mvrnorm_cpp(theta_current,sigma_proposal);
    
    //evaluate posterior of proposed theta
    double target_theta_propose;
    target_theta_propose = as<double>(wrap(target(theta_propose)));
    bool is_accepted;
    double log_acceptance;

    //if posterior is 0 then do not compute anything else and immediately reject
    if(!std::isfinite(target_theta_propose)){
      log_acceptance = -std::numeric_limits<double>::infinity();
    } else {
      //compute Metropolis-Hastings ratio (acceptance probability)
      log_acceptance = target_theta_propose - target_theta_current;
      // double A = as<double>(runif(0,1));
      double A = R::runif(0,1);
      if(log(A) < log_acceptance){ //accept the proposal
        is_accepted = true;
        theta_current = theta_propose;
      } else { //reject the proposal
        is_accepted = false;
      }
    }

    //update acceptance rate
    if(i==0){
      acceptance_rate = is_accepted;
    } else {
      acceptance_rate = acceptance_rate + (is_accepted - acceptance_rate) / i;
    }
    
    //update empirical covariance matrix 
    arma::vec residual = theta_current - theta_mean;
    theta_mean = theta_mean + (residual/i); //update empirical mean theta
    sigma_empirical = update_sigma(sigma_empirical,residual,i); //update empirical covariance matrix
    
    //store trace information
    theta_samp.row(i) = theta_current.t();
    sigma_trace.slice(i) = sigma_proposal;
    
  }
  
  //return output
  return(List::create(Named("theta_trace")=theta_samp,Named("sigma_trace")=sigma_trace,Named("acceptance_rate")=acceptance_rate));
}


