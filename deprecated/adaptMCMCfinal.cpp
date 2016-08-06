#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/*
 * Hi Yi-han, here's the MCMC code!
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
 * adaptMCMC returns an R list:
 * theta_trace a matrix recording the trace of the sampled parameters
 * sigma_empirical is an array recording the trace of the estimated covariance matrix of the target
 * acceptance_rate is a numeric that is the acceptance rate at the end of the MCMC run
 */
//[[Rcpp::export]]
List adaptMCMC(Function target, arma::vec init_theta, arma::mat covmat, int n_iterations, int adapt_size_start, int adapt_shape_start, int info, double adapt_size_cooling = 0.99){
  
  arma::vec theta_current = init_theta;
  arma::vec theta_propose = init_theta;
  arma::mat covmat_proposal = covmat;
  
  arma::mat covmat_proposal_init = covmat_proposal;
  bool adapting_size = false;
  bool adapting_shape = false;
  
  double target_theta_current = as<double>(wrap(target(theta_current)));

  arma::mat theta_trace = arma::zeros(n_iterations,init_theta.n_elem);
  arma::cube sigma_empirical = arma::zeros(covmat.n_rows,covmat.n_cols,n_iterations);

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
    
  }
  
  return(List::create(Named("theta_trace")=theta_trace,Named("sigma_empirical")=sigma_empirical,Named("acceptance_rate")=acceptance_rate));
}


/***R
#example function to sample from (https://en.wikipedia.org/wiki/Rosenbrock_function)
p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

#example of non-adaptive random walk MCMC in 2 dimensions
set.seed(123)
rw_output <- adaptMCMC(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),
                                n_iterations=1e3,
                                adapt_size_start=0,adapt_shape_start=0,info=1)

par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(rw_output$theta_trace, type='l')

matplot(rw_output$theta_trace,type="l")

par(mfrow=c(1,1))

#example of adaptive MCMC in 2 dimensions
set.seed(123)
adapt_output <- adaptMCMC(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),
                       n_iterations=1e3,
                       adapt_size_start=20,adapt_shape_start=50,info=1)

par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(adapt_output$theta_trace, type='l')

matplot(adapt_output$theta_trace,type="l")

par(mfrow=c(1,1))
*/