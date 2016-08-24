#include <Rcpp.h>
using namespace Rcpp;

/*
 * This file contains code relevant to a C++ implementation of a Gibbs sampler (a specific case of Metropolis-Hastings MCMC) to
 * sample from a truncated multivariate Gaussian distribution.
 */


// CDF_norm returns the cumulative distribution function of the standard normal (N(0,1)) evaluated at x (value)
// [[Rcpp::export]]
double CDF_norm(double x){
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

// invCDF_norm returns the inverse cumulative distribution function of the standard normal (N(0,1)) evaluated at x (probability)
// double invCDF_norm()
// 
// 
// double AffxStatistics::norminv(double p)
// {
//   const double a[6] = {
//     -3.969683028665376e+01,  2.209460984245205e+02,
//     -2.759285104469687e+02,  1.383577518672690e+02,
//     -3.066479806614716e+01,  2.506628277459239e+00
//   };
//   const double b[5] = {
//     -5.447609879822406e+01,  1.615858368580409e+02,
//     -1.556989798598866e+02,  6.680131188771972e+01,
//     -1.328068155288572e+01
//   };
//   const double c[6] = {
//     -7.784894002430293e-03, -3.223964580411365e-01,
//     -2.400758277161838e+00, -2.549732539343734e+00,
//     4.374664141464968e+00,  2.938163982698783e+00
//   };
//   const double d[4] = {
//     7.784695709041462e-03,  3.224671290700398e-01,
//     2.445134137142996e+00,  3.754408661907416e+00
//   };
//   
//   register double q, t, u;
//   
//   if ((p != p) || p > 1.0 || p < 0.0)
//     return NaN;
//   if (p == 0.0)
//     return NaN;
//   if (p == 1.0)
//     return  NaN;
//   q = min(p,1-p);
//   if (q > 0.02425) {
//     /* Rational approximation for central region. */
//     u = q-0.5;
//     t = u*u;
//     u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
//       /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
//   } else {
//     /* Rational approximation for tail region. */
//     t = sqrt(-2*log(q));
//     u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
//       /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
//   }
//   /* The relative error of the approximation has absolute value less
//    than 1.15e-9.  One iteration of Halley's rational method (third
//    order) gives full machine precision... */
//   t = normcdf(u)-q;    /* error */
//   t = t*sqrt(2*PI)*exp(u*u/2);   /* f(u)/df(u) */
//   u = u-t/(1+u*t/2);     /* Halley's method */
//   
//   return (p > 0.5 ? -u : u);
// };