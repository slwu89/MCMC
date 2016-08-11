library(tmvtnorm)

update_sigmaR <- function(cov_mat,residual,i){
  return((cov_mat * (i-1) + (i-1) / i * residual %*% t(residual)) / i)
}


adaptMCMC_R <- function (target, init.theta, n.iterations,
                         covmat = NULL,adapt.size.start = NULL, 
                         adapt.size.cooling = 0.99, adapt.shape.start = NULL,
                         info)
{
  
  theta.current <- init.theta
  theta.propose <- init.theta
  covmat.proposal <- covmat
  
  #store empirical covariance matrix and other information
  covmat.empirical.trace <- array(0,c(nrow(covmat),ncol(covmat),n.iterations))
  covmat.proposal.trace <- array(0,c(nrow(covmat),ncol(covmat),n.iterations))
  residual.trace <- matrix(0,nrow=n.iterations,ncol=length(init.theta))
  theta.mean.trace <- matrix(0,nrow=n.iterations,ncol=length(init.theta))
  scaling.sd.trace <- rep(0,n.iterations)
  acceptance.trace <- rep(0,n.iterations)
  target.theta.current.trace <- rep(0,n.iterations)
  
  covmat.proposal.init <- covmat.proposal
  adapting.size <- FALSE
  adapting.shape <- FALSE
  
  target.theta.current <- target(theta.current)
  
  trace <- data.frame(t(theta.current))
  acceptance.rate <- 0
  scaling.sd <- 1
  covmat.empirical <- covmat.proposal
  covmat.empirical[, ] <- 0
  theta.mean <- theta.current
  
  for (i.iteration in seq_len(n.iterations)) {
    
    message("new iteration")
    
    if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start &&
        acceptance.rate * i.iteration < adapt.shape.start) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      }
      scaling.sd <- scaling.sd * exp(adapt.size.cooling^(i.iteration - adapt.size.start) * (acceptance.rate - 0.234))
      print(paste("scaling.sd:",scaling.sd,"covmat.proposal:"))
      covmat.proposal <- scaling.sd^2 * covmat.proposal.init
      print(covmat.proposal)
    } else if (!is.null(adapt.shape.start) && acceptance.rate *
               i.iteration >= adapt.shape.start) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        adapting.shape <- TRUE
      }
      covmat.proposal <- 2.38^2/length(init.theta) * covmat.empirical
      print(paste("covmat.proposal:"))
      print(covmat.proposal)
    }
    
    if (any(diag(covmat.proposal) < .Machine$double.eps)) {
      print(covmat.proposal)
      stop("non-positive definite covmat", call. = FALSE)
    }
    
    theta.propose <- as.vector(rmvnorm(1,mean = theta.current, sigma = covmat.proposal))
    target.theta.propose <- target(theta.propose)
    print(paste("theta.propose and target.theta.proposal:"))
    print(theta.propose)
    print(target.theta.propose)
    
    
    if (!is.finite(target.theta.propose)) {
      log.acceptance <- -Inf
    }
    else {
      log.acceptance <- target.theta.propose - target.theta.current
      print(paste("log.acceptance step 1",log.acceptance))
      log.acceptance <- log.acceptance + dmvnorm(x = theta.current,mean = theta.propose,sigma = covmat.proposal,log = TRUE)
      print(paste("log.acceptance step 2",log.acceptance))
      log.acceptance <- log.acceptance - dmvnorm(x = theta.propose,mean = theta.current,sigma = covmat.proposal,log = TRUE)
      print(paste("log.acceptance step 3",log.acceptance))
    }
    if (log(runif(1)) < log.acceptance) {
      theta.current <- theta.propose
      target.theta.current <- target.theta.propose
      is.accepted <- TRUE
    } else {
      is.accepted <- FALSE
    }
    trace <- rbind(trace, c(theta.current))
    acceptance.rate <- acceptance.rate + (is.accepted - acceptance.rate)/i.iteration
    
    residual <- as.vector(theta.current-theta.mean)
    print(paste("residual:"))
    print(residual)
    covmat.empirical <- (covmat.empirical*(i.iteration-1)+(i.iteration-1)/i.iteration*residual%*%t(residual))/i.iteration
    print(paste("covmat.empirical"))
    print(covmat.empirical)
    theta.mean <- theta.mean + residual/i.iteration
    print(paste("theta.mean:",theta.mean))
    
    #store empirical covariance matrix
    covmat.empirical.trace[,,i.iteration] <- covmat.empirical
    residual.trace[i.iteration,] <- residual
    theta.mean.trace[i.iteration,] <- theta.mean
    scaling.sd.trace[i.iteration] <- scaling.sd
    target.theta.current.trace[i.iteration] <- target.theta.current
    covmat.proposal.trace[,,i.iteration] <- covmat.proposal
    acceptance.trace[i.iteration] <- acceptance.rate
    scaling.sd.trace[i.iteration] <- scaling.sd
    
    if(i.iteration %% info == 0){
      print(paste0("Chain information at iter: ",i.iteration," acceptance.rate: ",acceptance.rate))
    }
    
  }
  return(list(theta_trace = trace, acceptance.rate = acceptance.rate,
              covmat.empirical.trace = covmat.empirical.trace,residual.trace=residual.trace,
              theta.mean.trace=theta.mean.trace,
              covmat.proposal.trace=covmat.proposal.trace,target.theta.current.trace=target.theta.current.trace,
              acceptance.trace=acceptance.trace,scaling.sd.trace=scaling.sd.trace))
}

p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

par(mfrow=c(1,2))

#test random walk mcmc
set.seed(123)
banana_out <- adaptMCMC_R(target=p.log,init.theta=c(10,10),covmat=diag(c(1,1)),
                          adapt.size.start=10,adapt.shape.start=20,
                          n.iterations=1e3,info=1)

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_out$theta_trace, type='l')

matplot(banana_out$theta_trace,type="l")

par(mfrow=c(1,1))

