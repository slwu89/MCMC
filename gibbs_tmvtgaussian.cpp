#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/*
 * This file contains code relevant to a C++ implementation of a Gibbs sampler (a specific case of Metropolis-Hastings MCMC) to
 * sample from a truncated multivariate Gaussian distribution.
 */

// sub1 returns a matrix x[-i,-i]
// [[Rcpp::export]]
arma::mat sub1(arma::mat x, int i) {
  x.shed_col(i);
  x.shed_row(i);
  return x;
}

// sub2 returns a matrix x[a,-b]
// [[Rcpp::export]]
arma::mat sub2(arma::mat x, int a, int b){
  x.shed_col(b);
  return(x.row(a));
}

// negSubCol returns a column vector x[-i]
// [[Rcpp::export]]
arma::vec negSubCol(arma::vec x, int i){
  x.shed_row(i);
  return(x);
}

// negSubRow returns a row vector x[-i]
// [[Rcpp::export]]
arma::rowvec negSubRow(arma::rowvec x, int i){
  x.shed_col(i);
  return(x);
}

// r8poly_value_horner evaluates a polynomial using Horner's method
// [[Rcpp::export]]
double r8poly_value_horner(int m, double c[], double x){
  
  double value;
  value = c[m];
 
  for(int i = m - 1; 0 <= i; i--){
    value = value * x + c[i];
  }
  
  return value;
}


/*
 * CDF_norm returns the cumulative distribution function of the standard normal N(mu,sigma) evaluated at x (value);
 * this function replicates the functionality of pnorm from base R
 */
// [[Rcpp::export]]
double CDF_norm(double x_input, double mu, double sigma){
  
  //map input to N(0,1)
  double x = (x_input - mu) / sigma;
  
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;

  if (fabs(x) <= 1.28){
    y = 0.5 * x * x;
    q = 0.5 - fabs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))));
  } else if (fabs(x) <= 12.7){
    y = 0.5 * x * x;
    q = exp (-y) * b0 / (fabs(x) - b1 + b2  / (fabs(x) + b3 + b4  / (fabs(x) - b5 + b6  / (fabs (x) + b7 - b8  / (fabs(x) + b9 + b10 / (fabs(x) + b11))))));
  } else {
    q = 0.0;
  }

  if (x < 0.0)
  {
    cdf = q;
  } else {
    cdf = 1.0 - q;
  }
  
  return cdf;
}


/*
 * invCDF_norm returns the inverse cumulative distribution function of the standard normal N(mu,sigma) evaluated at x (probability);
 * this function replicates the functionality of qnorm from base R
 */
// [[Rcpp::export]]
double invCDF_norm(double p, double mu, double sigma){
  
  double out;
  
  //map input to N(0,1)
  if(p < 0.0 || 1.0 < p){
    Rcout << "\n" << std::endl;
    Rcout << "NORMAL_MS_CDF_INV - Fatal error!\n" << std::endl;
    Rcout << "  CDF < 0 or 1 < CDF.\n" << std::endl;
    return(0.0);
  }
  
  double a[8] = {
    3.3871328727963666080,     1.3314166789178437745E+2,
    1.9715909503065514427E+3,  1.3731693765509461125E+4,
    4.5921953931549871457E+4,  6.7265770927008700853E+4,
    3.3430575583588128105E+4,  2.5090809287301226727E+3 };
  double b[8] = {
    1.0,                       4.2313330701600911252E+1,
    6.8718700749205790830E+2,  5.3941960214247511077E+3,
    2.1213794301586595867E+4,  3.9307895800092710610E+4,
    2.8729085735721942674E+4,  5.2264952788528545610E+3 };
  double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770E-1,
    2.27238449892691845833E-2,  7.74545014278341407640E-4 };
  double const1 = 0.180625;
  double const2 = 1.6;
  double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550E-1,
    1.48103976427480074590E-1,  1.51986665636164571966E-2,
    5.47593808499534494600E-4,  1.05075007164441684324E-9 };
  double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230E-1,
    2.65321895265761230930E-2,  1.24266094738807843860E-3,
    2.71155556874348757815E-5,  2.01033439929228813265E-7 };
  double f[8] = {
    1.0,                        5.99832206555887937690E-1,
    1.36929880922735805310E-1,  1.48753612908506148525E-2,
    7.86869131145613259100E-4,  1.84631831751005468180E-5,
    1.42151175831644588870E-7,  2.04426310338993978564E-15 };
  double q;
  double r;
  double split1 = 0.425;
  double split2 = 5.0;
  double value;
  
  if (p <= 0.0){
    value = -HUGE_VAL;
    return value;
  }
  if (1.0 <= p){
    value = HUGE_VAL;
    return value;
  }

  q = p - 0.5;
  
  if(fabs(q) <= split1){
    r = const1 - q * q;
    value = q * r8poly_value_horner(7,a,r) / r8poly_value_horner(7,b,r);
  } else {
    if(q < 0.0){
      r = p;
    } else {
      r = 1.0 - p;
    }
     if (r <= 0.0){
      value = HUGE_VAL;
    } else {
      r = sqrt ( - log ( r ) );
      if(r <= split2){
        r = r - const2;
        value = r8poly_value_horner(7,c,r) / r8poly_value_horner(7,d,r);
      } else {
        r = r - split2;
        value = r8poly_value_horner(7,e,r) / r8poly_value_horner(7,f,r);
      }
    }
    
    if (q < 0.0){
      value = - value;
    }
  }
  
  out = mu + sigma * value;
  return(out);
}


/*
 * rtnorm_gibbs returns a sample of size n from the specificed truncated Gaussian distribution.
 * n: integer number of samples to take
 * mu: mean of distribution
 * sigma: standard deviation of distribution
 * a: lower truncation bound
 * b: upper truncation bound
 */
// [[Rcpp::export]]
NumericVector rtnorm_gibbs(int n, double mu, double sigma, double a, double b){

  //sample from uniform distribution on unit interval
  NumericVector F = runif(n);

  //Phi(a) and Phi(b)
  double Fa = CDF_norm(a,mu,sigma);
  double Fb = CDF_norm(b,mu,sigma);
  
  NumericVector F_out(F.length());
  for(int i=0; i < F.length(); i++){
    double p_i = F[i] * (Fb - Fa) + Fa;
    F_out[i] = invCDF_norm(p_i,mu=0.0,sigma=1.0);
  }
  
  NumericVector out(F.length());
  for(int i=0; i < out.length(); i++){
    out[i] = mu + sigma * F_out[i];
  }

  return(out);
}

// [[Rcpp::export]]
void test_gibbs(){
  NumericVector gibbs = rtnorm_gibbs(5,0.0,1.0,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());
  Rcout << gibbs << std::endl;
  
}

/***R
library(tmvtnorm)
set.seed(123)
tmvtnorm:::rtnorm.gibbs(n=10)
set.seed(123)
rtnorm_gibbs(10,0.0,1.0,-Inf,Inf)
set.seed(123)
test_gibbs()
*/


// /*
//  * rtmvnorm_gibbs returns a sample of size n from the specified truncated multivariate Gaussian distribution
//  * 
//  */
// arma::mat rtmvnorm_gibbs(int n, arma::vec mu, arma::mat sigma, arma::vec lower, arma::vec upper, arma::vec init_state){
//   
//   int d = mu.n_elem; //check dimension of target distribution
//   arma::mat trace = arma::zeros(n,d); //trace of MCMC chain
//   
//   //draw from U(0,1)
//   NumericVector U = runif(n*d);
//   int l = 0; //iterator for U
//   
//   //calculate conditional standard deviations
//   arma::vec sd(d);
//   arma::cube P = arma::zeros(1,d-1,d);
//   
//   for(int i=0; i<d; i++){
//     //partitioning of sigma
//     arma::mat Sigma = sub1(sigma,i);
//     double sigma_ii = sigma(i,i);
//     arma::rowvec Sigma_i = sub2(sigma,i,i);
// 
//     P.slice(i) = Sigma_i * Sigma.i();
//     double p_i = Rcpp::as<double>(wrap(P.slice(i) * Sigma_i.t()));
//     sd(i) = sqrt(sigma_ii - p_i);
//   }
//   
//   arma::vec x = init_state;
//   
//   //run Gibbs sampler for specified chain length (MCMC chain of n samples)
//   for(int j=0; j<n; j++){
//     
//     //sample all conditional distributions
//     for(int i=0; i<d; i++){
//       
//       //calculation of conditional expectation and conditional variance
//       arma::rowvec slice_i = P.slice(i);
//       arma::vec slice_i_times = slice_i * (negSubCol(x,i) - negSubCol(x,i));
//       double slice_i_times_double = Rcpp::as<double>(wrap(slice_i_times));
//       double mu_i = mu(i) + slice_i_times_double;
//       
//       //transformation
//       double Fa = CDF_norm(lower(i),mu_i,sd(i));
//       double Fb = CDF_norm(upper(i),mu_i,sd(i));
//       x(i) = mu_i + sd(i) * invCDF_norm(U(l) * (Fb - Fa) + Fa,mu=0.0,sigma=1.0);
//       l = l + 1;
//       
//     }
//     
//     trace.row(j) = x.t();
//     
//   }
//   
//   return(trace);
// }

