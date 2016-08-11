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
  for(int i=1; i<=iterations; i++){  //REMEMBER TO SUBSET BY i-1 !!!!!!!!!!!!!!!
    
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
  
  
/***R
memory.limit(size=540000)

p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

set.seed(123)
mcmc1 <- adaptMCMC(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),iterations=1e3,
                         adapt_size_start=10,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=20,
                         info=1e2)

par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(mcmc1$theta_trace, type='l')

matplot(mcmc1$theta_trace,type="l")

par(mfrow=c(1,1))


#############################################
####plot the empirical estimation of sigma###
#############################################

# load packages
# library(foreach)
# library(doSNOW)
# library(parallel)

#estimate kernel denstiy of empirical sigma
# cl <- makeCluster(spec=detectCores()-2)
# registerDoSNOW(cl)
# sigma_kdens <- foreach(i=1:nrow(mcmc1$theta_trace),.packages=c("MASS"),.verbose=TRUE) %dopar% {
#  if(sum(mcmc1$sigma_trace[,,i])==0){
#    return(NULL)
#  }
#  reps <- mvrnorm(1e4,mu=mcmc1$theta_trace[i,],Sigma=mcmc1$sigma_trace[,,i])
#  dd <- kde2d(reps[,1],reps[,2],n=200)
#  return(dd)
# }
# 
# stopCluster(cl)
# rm(cl)

# sigma_kdens <- sigma_kdens[!sapply(sigma_kdens,is.null)]

#generate plots
# bounds <- sapply(sigma_kdens,function(x){
#  minX=min(x$x);
#  maxX=max(x$x);
#  minY=min(x$y);
#  maxY=max(x$y);
#  return(c(minX,maxX,minY,maxY))
# })
# 
# minX <- min(t(bounds)[,1])
# maxX <- max(t(bounds)[,2])
# minY <- min(t(bounds)[,3])
# maxY <- max(t(bounds)[,4])

#need to figure out how to "fill in" the rest of the plotting area with red
# for(i in 1:length(sigma_kdens)){
#   image(sigma_kdens[[i]],xlim=c(minX,maxX),ylim=c(minY,maxY),useRaster=T)
#   contour(x=sigma_kdens[[i]]$x,y=sigma_kdens[[i]]$y,z=sigma_kdens[[i]]$z,add=T,lty=2)
# }

# file_root <- "C:/Users/Administrator/Dropbox/GitHub/MCMC/graphics/"
# 
# for(i in 1:length(sigma_kdens)){
#  file_path <- file.path(paste0(file_root,"sigma_plot",i,".jpg"))
#  jpeg(file_path,quality=100,width=640,height=640)
#  image(sigma_kdens[[i]],useRaster=T)
#  contour(x=sigma_kdens[[i]]$x,y=sigma_kdens[[i]]$y,z=sigma_kdens[[i]]$z,add=T,lty=2)
#  dev.off()
# }
 
# ggplot(data=melt(sigma_kdens[[1]]$z),aes(x=Var1,y=Var2,fill=value)) +
#   geom_raster() +
#   geom_contour(aes(z=value)) +
#   guides(fill=FALSE) +
#   theme_bw()
*/
  
  
  
/*
 * Simple Random Walk MCMC with Adaptive Proposal Kernel 
 * adaptMCMC returns an R list:
 * theta_trace a matrix recording the trace of the sampled parameters
 * sigma_empirical is an array recording the trace of the estimated covariance matrix of the target
 * acceptance_rate is a numeric that is the acceptance rate at the end of the MCMC run
 */
//[[Rcpp::export]]
List adaptMCMC_simple(Function target, arma::vec init_theta, arma::mat covmat, int n_iterations, int adapt_size_start, int adapt_shape_start, int info, double adapt_size_cooling = 0.99){
  
  arma::vec theta_current = init_theta;
  arma::vec theta_propose = init_theta;
  arma::mat covmat_proposal = covmat;
  
  arma::mat covmat_proposal_init = covmat_proposal;
  bool adapting_size = false;
  bool adapting_shape = false;
  
  double target_theta_current = as<double>(wrap(target(theta_current)));
  
  //henceforth 
  //target.theta.current$log.density <- target.theta.current
  //target.theta.current$trace <- theta.current
  
  arma::mat theta_trace = arma::zeros(n_iterations,init_theta.n_elem);
  arma::cube sigma_empirical = arma::zeros(covmat.n_rows,covmat.n_cols,n_iterations);
  
  //store extra stuff for debugging
  arma::vec acceptance_trace = arma::zeros(n_iterations);
  arma::vec target_theta_current_trace = arma::zeros(n_iterations);
  arma::vec scaling_sd_trace = arma::zeros(n_iterations);
  arma::cube covmat_proposal_trace = arma::zeros(covmat.n_rows,covmat.n_cols,n_iterations);
  arma::mat residual_trace = arma::zeros(n_iterations,init_theta.n_elem);
  arma::mat theta_mean_trace = arma::zeros(n_iterations,init_theta.n_elem);
  
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
    acceptance_trace(i-1) = acceptance_rate;
    scaling_sd_trace(i-1) = scaling_sd;
    target_theta_current_trace(i-1) = target_theta_current;
    covmat_proposal_trace.slice(i-1) = covmat_proposal;
    residual_trace.row(i-1) = residual.t();
    theta_mean_trace.row(i-1) = theta_mean.t();
  }
  
  // return(List::create(Named("theta_trace")=theta_trace,Named("sigma_empirical")=sigma_empirical,Named("acceptance_rate")=acceptance_rate));
  return(List::create(Named("theta_trace")=theta_trace,Named("sigma_empirical")=sigma_empirical,Named("acceptance_rate")=acceptance_rate,Named("acceptance_trace")=acceptance_trace,Named("target_theta_current_trace")=target_theta_current_trace,Named("covmat_proposal_trace")=covmat_proposal_trace,Named("residual_trace")=residual_trace,Named("theta_mean_trace")=theta_mean_trace,Named("scaling_sd_trace")=scaling_sd_trace));
}

/***R
p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

set.seed(123)
mcmc2 <- adaptMCMC_simple(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),
                                n_iterations=1e3,
                         adapt_size_start=10,adapt_shape_start=20,info=1)
  
  par(mfrow=c(1,2))
  
  x1 <- seq(-15, 15, length=100)
  x2 <- seq(-15, 15, length=100)
  d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
  image(x1, x2, exp(d.banana), col=cm.colors(60))
  contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
  lines(mcmc2$theta_trace, type='l')
  
  matplot(mcmc2$theta_trace,type="l")
  
  par(mfrow=c(1,1))
  
*/
  
  
