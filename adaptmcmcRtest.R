update_sigmaR <- function(cov_mat,residual,i){
  return((cov_mat * (i-1) + (i-1) / i * residual %*% t(residual)) / i)
}

adaptMCMC <- function(target,theta_init,sigma,cooling,adapt_size,adapt_shape,iterations,info,max_scaling_sd=50){
  
  theta_current = theta_init
  theta_propose = theta_init
  
  sigma_proposal = sigma
  sigma_proposal_init = sigma
  
  adapting_size <- FALSE
  adapting_shape <- FALSE
  
  target_theta_current <- target(theta_current)
  
  theta_samp <- matrix(data=NA,nrow=iterations,ncol=length(theta_init))
  sigma_trace <- array(NA,dim=c(length(theta_init),length(theta_init),iterations))
  
  acceptance_rate = 0.0
  scaling_sd = 1.0
  scaling_multiplier = 1.0
  
  sigma_empirical <- matrix(0,nrow(sigma),ncol(sigma))
  theta_mean <- theta_current
  
  for(i in 1:iterations){
    
    if(i >= adapt_size & acceptance_rate*i < adapt_shape){
      if(!adapting_size){
        print(paste("Begin adapting size of sigma at iteration",i))
        adapting_size <- TRUE
      }
      scaling_multiplier <- exp(cooling^(i-adapt_size) * (acceptance_rate - 0.234))
      scaling_sd <- scaling_sd * scaling_multiplier
      scaling_sd <- min(c(scaling_sd,max_scaling_sd))
      sigma_proposal_new <- scaling_sd^2 * sigma_proposal_init
      if(!any(diag(sigma_proposal_new) < 2e-16)){
        sigma_proposal <- sigma_proposal_new
      }
    } else if(acceptance_rate*i >= adapt_shape){
      if(!adapting_shape){
        print(paste("begin adapting shape at",i))
        adapting_shape <- TRUE
      }
      scaling_sd <- 2.38/sqrt(length(theta_init))
      sigma_proposal <- scaling_sd^2 * sigma_empirical
    }
    
    if(i %% info == 0){
      print(paste("at iter",i,"acc rate",acceptance_rate))
    }
    
    theta_propose = mvrnorm_samp(theta_current,sigma_proposal)
    
    target_theta_propose = target(theta_propose)
    
    if(!is.finite(target_theta_propose)){
      log_acceptance <- -Inf
    } else {
      log_acceptance <- target_theta_propose - target_theta_current
      log_acceptance = log_acceptance + mvrnorm_pdf(theta_current,theta_propose,sigma_proposal);
      log_acceptance = log_acceptance - mvrnorm_pdf(theta_propose,theta_current,sigma_proposal);
      if(log(runif(1)) < log_acceptance){
        is_accepted <- TRUE
        theta_current <- theta_propose
        target_theta_current <- target_theta_propose
      } else {
        is_accepted <- FALSE
      }
    }
    
    if(i==1){
      acceptance_rate = is_accepted;
    } else {
      acceptance_rate = acceptance_rate + (is_accepted - acceptance_rate) / i;
    }
    
    residual = theta_current - theta_mean
    print(paste("residual",residual))
    sigma_empirical = update_sigmaR(sigma_empirical,residual,i)
    print(paste("sigma_empirical",sigma_empirical))
    theta_mean = theta_mean + residual/i
    print(paste("theta_mean",theta_mean))
    
    theta_samp[i,] <- theta_current
    sigma_trace[,,i] <- sigma_proposal
    

  }
  return(list("theta_trace"=theta_samp,"sigma_trace"=sigma_trace))
}
