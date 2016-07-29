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
set.seed(123)
banana_out <- rw_mcmc(target=p.log,theta_init=c(10,10),sigma=diag(c(1,1)),iterations=1e3,info=1e2)

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_out$trace, type='l')

matplot(banana_out$trace,type="l")

#test adaptive mcmc
set.seed(123)
banana_adapt_out <- adapt_mcmc(target=p.log,theta_init=c(10,10),sigma=diag(c(1,1)),cooling=0.99,adapt_size=10,adapt_shape=20,iterations=1e3,info=1e2)

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_adapt_out$theta_trace, type='l')

matplot(banana_adapt_out$theta_trace,type="l")

par(mfrow=c(1,1))


#try fitR package version of adaptive mcmc OR adaptmcmcRtest.R version
library(fitR)
p.log.fitR <- function(x) {
  if(any(names(x) %in% c("par1","par2"))){
    x1 <- x[["par1"]]
    x2 <- x[["par2"]]
  } else {
    x1 <- x[1]
    x2 <- x[2]
  }
  B <- 0.03 # controls 'bananacity'
  -x1^2/200 - 1/2*(x2+B*x1^2-100*B)^2
}

covmat_fitR <- matrix(c(.00001,0,0,10),2,2,dimnames=list(c("par1","par2"),c("par1","par2")))

set.seed(123)
banana_fitR_out <- mcmcMH(target=p.log.fitR,init.theta=c(par1=10,par2=10),covmat=covmat_fitR,adapt.size.start=10,adapt.shape.start=20,n.iterations=1e3)
banana_fitR_out <- adaptMCMC(target=p.log,theta_init=c(10,10),sigma=diag(c(1,1)),cooling=0.99,adapt_size=10,adapt_shape=20,iterations=1e3,info=1e2)


par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log.fitR), nrow=100)
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


###testing Adaptive MCMC only###

#testing with Goldenstein-Price function
goldpr <- function(xx){
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1a <- (x1 + x2 + 1)^2
  fact1b <- 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2
  fact1 <- 1 + fact1a*fact1b
  
  fact2a <- (2*x1 - 3*x2)^2
  fact2b <- 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2
  fact2 <- 30 + fact2a*fact2b
  
  y <- fact1*fact2
  return(-log(y))
}

goldpr_adapt_out <- adapt_mcmc(target=goldpr,theta_init=c(1,1),sigma=diag(c(1,1)),cooling=0.99,adapt_size=1e2,adapt_shape=1.1e2,iterations=1e3,info=1e2)

par(mfrow=c(1,2))

x1 <- seq(-2, 2, length=100)
x2 <- seq(-2, 2, length=100)
d.goldpr <- matrix(apply(expand.grid(x1, x2), 1, goldpr), nrow=100)
image(x1, x2, d.goldpr, col=cm.colors(100))
contour(x1, x2, d.goldpr, add=TRUE, col=gray(0.6))
lines(goldpr_adapt_out$theta_trace, type='l')

matplot(goldpr_adapt_out$theta_trace,type="l")

par(mfrow=c(1,1))

#testing with the Banana function
banana <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

banana_adapt_out <- adapt_mcmc(target=p.log,theta_init=c(10,10),sigma=diag(c(1,1)),cooling=0.5,adapt_size=5e2,adapt_shape=1e5,iterations=5e3,info=1e2)

par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_adapt_out$theta_trace, type='l')

matplot(banana_adapt_out$theta_trace,type="l")

par(mfrow=c(1,1))


###testing the sigma_empirical updating function
cov_mat <- matrix(c(5,0,0,0.1),2,2)
residual <- c(.25,1.5)
i <- 5

update_sigmaR <- function(cov_mat,residual,i){
  out <- (cov_mat*(i-1)+(i-1)/i*residual%*%t(residual))/i
  return(out)
}

update_sigmaR(cov_mat,residual,i)

(cov_mat*(i-1)+(i-1)/i*residual%*%t(residual))/i

cppFunction(code="arma::mat update_sigmaC(arma::mat cov_mat, arma::vec residual, int i){
  arma::mat out = (cov_mat * (i-1) + (i-1)/i*residual * trans(residual))/i;
  return(out);
}",depends=c("RcppArmadillo"))

update_sigmaR(cov_mat,residual,i)
update_sigmaC(cov_mat,residual,i)


cov_mat <- matrix(rnorm(4),2,2)
theta_mean <- c(9.2778,7.9439)
theta_current <- c(10.6343,8.7317)
i <- 10

residual <- as.vector(theta_current-theta_mean)
(cov_mat*(i-1)+(i-1)/i*residual%*%t(residual))/i

covUpdateTest(cov_mat,residual,i)
update_sigma(cov_mat,residual,i)
