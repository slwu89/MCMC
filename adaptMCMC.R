############################################################
######Markov Chain Monte Carlo in R and C++ (via Rcpp)######
############################################################

library(Rcpp)
library(RcppArmadillo)

Rcpp::sourceCpp('C:/Users/WuS/Dropbox/GitHub/R_repo/adaptMCMC.cpp')

#function to test on 
p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

par(mfrow=c(2,2))

#test random walk mcmc
banana_out <- rw_mcmc(target=p.log,theta_init=c(10,10),sigma=diag(c(1,1)),iterations=2e4,info=1e3)

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_out$trace, type='l')

matplot(banana_out$trace,type="l")

#test adaptive mcmc
banana_adapt_out <- adapt_mcmc(target=p.log,theta_init=c(10,10),sigma=diag(c(1,1)),cooling=0.99,adapt_size=10,adapt_shape=20,iterations=5e5,info=1e4)

image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_adapt_out$theta_trace, type='l')

matplot(banana_adapt_out$theta_trace,type="l")

par(mfrow=c(1,1))


#try fitR package version of adaptive mcmc
p.log.fitR <- function(x) {
  x1 <- x[["par1"]]
  x2 <- x[["par2"]]
  B <- 0.03 # controls 'bananacity'
  -x1^2/200 - 1/2*(x2+B*x1^2-100*B)^2
}

covmat_fitR <- matrix(c(1,0,0,1),2,2,dimnames=list(c("par1","par2"),c("par1","par2")))

banana_fitR_out <- mcmcMH(target=p.log.fitR,init.theta=c(par1=10,par2=10),covmat=covmat_fitR,adapt.size.start=250,adapt.shape.start=500,n.iterations=2e4)

par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_fitR_out$trace, type='l')

matplot(banana_fitR_out$trace,type="l")

par(mfrow=c(1,1))

##sink the function output of adaptive mcmc for debugging

# sink("debug.txt")
# banana_adapt_out <- adapt_mcmc(target=p.log,theta_init=c(-10,10),sigma=diag(c(5,0.1)),cooling=0.99,adapt_size=250,adapt_shape=500,iterations=3e3,info=100)
# unlink("debug.txt")
# closeAllConnections()



###testing the sigma_empirical updating function
cov_mat <- matrix(c(5,0,0,0.1),2,2)
residual <- c(0,0)
i <- 2

update_sigmaR <- function(cov_mat,residual,i){
  out <- (cov_mat*(i-1)+(i-1)/i*residual%*%t(residual))/i
  return(out)
}

update_sigmaR(cov_mat,residual,i)

update_sigma(cov_mat,residual,i)

cppFunction(code="arma::mat update_sigmaC(arma::mat cov_mat, arma::vec residual, int i){
  arma::mat out = (cov_mat * (i-1) + (i+1)/i*residual * residual.t())/i;
  return(out);
}",depends=c("RcppArmadillo"))

update_sigmaC(cov_mat,residual,i)
