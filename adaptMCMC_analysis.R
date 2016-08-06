##############################################
######Analysis of Adaptive MCMC routines######
######Sean Wu 8/5/2016########################
##############################################

Rcpp::sourceCpp('/Users/slwu89/Dropbox/GitHub/MCMC/adaptMCMC_source.cpp')


##################################
###Rosenbrock (banana) Function###
##################################

p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

#plot output from complex sampler
set.seed(123)
banana1 <- adaptMCMC(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),iterations=1e3,
                   adapt_size_start=10,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=20,
                   info=1e2)

par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana1$theta_trace, type='l')

matplot(banana1$theta_trace,type="l")

par(mfrow=c(1,1))

#plot output from simple sampler
set.seed(123)
banana2 <- adaptMCMC_simple(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),
                          n_iterations=1e3,
                          adapt_size_start=10,adapt_shape_start=20,verbose=FALSE,info=1e2)

par(mfrow=c(1,2))

image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana2$theta_trace, type='l')

matplot(banana2$theta_trace,type="l")

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


################################
###Goldenstein-Price Function###
################################

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

goldpr1 <- adaptMCMC(target=goldpr,init_theta=c(1,1),covmat=diag(c(1,1)),iterations=1e3,
             adapt_size_start=10,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=20,
             info=1e2)

par(mfrow=c(1,2))

x1 <- seq(-2, 2, length=100)
x2 <- seq(-2, 2, length=100)
d.goldpr <- matrix(apply(expand.grid(x1, x2), 1, goldpr), nrow=100)
image(x1, x2, d.goldpr, col=cm.colors(100))
contour(x1, x2, d.goldpr, add=TRUE, col=gray(0.6))
lines(goldpr1$theta_trace, type='l')

matplot(goldpr1$theta_trace,type="l")

par(mfrow=c(1,1))
