#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//This file include functions relevant to the truncated multivariate gaussian distribution


// mvrGaussian_density evaluates the PDF of the multivariate Gaussian distribution
// plase note this function returns the un-logged PDF to avoid infinite values
// [[Rcpp::export]]
double mvrGaussian_pdf(arma::colvec x, arma::colvec mu, arma::mat sigma){
  
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
  return(out);
}


// tmvrGaussian_density evaluates the PDF of the truncated multivariate Gaussian distribution with specified  bounds
// please note this function returns the un-logged PDF to avoid infinite values
// [[Rcpp::export]]
double tmvrGaussian_pdf(NumericVector x, NumericVector mu, NumericMatrix sigma, NumericVector lower, NumericVector upper){
  
  int n = x.size();
  double out;
  
  //generate boolean vector indicating if points in x fall inside support region
  LogicalVector support(n);
  for(int i=0; i < n; i++){
    if((x(i) >= lower[i]) && (x(i) <= upper[i])){
      support[i] = true;
    } else {
      support[i] = false;
    }
  }
  
  //check if any points fell outside the support region
  bool outside_support = false;
  for(int i=0; i < support.size(); i++){
    if(support[i] == false){
      outside_support = true;
    }
  }
  
  //if points fall outside support region immediately break and return 0
  if(outside_support){
    out = 0;
    return(out);
  }
  
  //if points fall inside support region calculate the density
  Environment mvtnorm("package:mvtnorm");
  Function pmvnorm = mvtnorm["pmvnorm"];
  double pdf;
  SEXP cdf_r;
  pdf = mvrGaussian_pdf(as<arma::colvec>(x),as<arma::colvec>(mu),as<arma::mat>(sigma));
  cdf_r = pmvnorm(lower,upper,mu,R_NilValue,sigma);
  double cdf = as<double>(cdf_r);
  out = pdf / cdf;
  
  return(out);
}


// // tmvrGaussian_gibbs samples from the truncated multivariate Gaussian distribution using Gibbs sampling
// // [[Rcpp::export]]
// arma::rowvec tmvrGaussian_gibbs(int n, NumericVector mu, NumericMatrix sigma, NumericVector lower, NumericVector upper, int burn_in, int thinning){
//   
// }

  
  
