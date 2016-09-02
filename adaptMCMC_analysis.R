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
banana1 <- adaptMCMC(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),n_iterations=1e3,
                   adapt_size_start=10,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=20,
                   verbose=FALSE,info=1e2)

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

#load packages
library(foreach)
library(doSNOW)
library(parallel)
library(ggplot2)
library(reshape2)
library(gganimate)

#estimate kernel denstiy of empirical sigma
cl <- makeCluster(spec=detectCores()-2)
registerDoSNOW(cl)
sigma_kdens <- foreach(i=1:nrow(banana1$theta_trace),.packages=c("MASS"),.verbose=TRUE) %dopar% {
 if(sum(banana1$sigma_empirical[,,i])==0){
   return(NULL)
 }
 reps <- mvrnorm(1e4,mu=banana1$theta_trace[i,],Sigma=banana1$sigma_empirical[,,i])
 dd <- kde2d(reps[,1],reps[,2],n=200)
 return(dd)
}

stopCluster(cl)
rm(cl)

sigma_kdens <- sigma_kdens[!sapply(sigma_kdens,is.null)]

#return melted data
get_animationData <- function(sigmaDat){
  
  cl <- makeCluster(spec=detectCores()-2)
  registerDoSNOW(cl)
  
  out <- foreach(i=1:length(sigmaDat),.packages=c("reshape2"),.verbose=TRUE,.combine="rbind") %dopar% {
    dat_i <- melt(sigmaDat[[i]]$z)
    dat_i$iter <- i
    return(dat_i)
  }
  
  stopCluster(cl)
  rm(cl)
  
  return(out)
}

###GGanimate Version###

sigma_aniDat <- get_animationData(sigma_kdens[10:15])

ggplot(data=sigma_aniDat[sigma_aniDat$iter==1,],aes(x=Var1,y=Var2,z=value,fill=value)) +
  geom_raster() +
  geom_contour() +
  guides(fill=FALSE) +
  theme_bw() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.title=element_blank(),title=element_blank(),panel.border=element_blank())

sigma_ggAni <- ggplot(data=sigma_aniDat,aes(x=Var1,y=Var2,z=value,fill=value,frame=iter)) +
  geom_raster() +
  geom_contour() +
  guides(fill=FALSE) +
  theme_bw() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.title=element_blank(),title=element_blank(),panel.border=element_blank())

gg_animate(sigma_ggAni)

###Base plotting with Animation package###

max_x <- max(sapply(sigma_kdens,function(x){max(x$x)}))
min_x <- min(sapply(sigma_kdens,function(x){min(x$x)}))
max_y <- max(sapply(sigma_kdens,function(x){max(x$y)}))
min_y <- min(sapply(sigma_kdens,function(x){min(x$y)}))

library(animation)

makeplot <- function(){
  for(i in 1:500){
    image(sigma_kdens[[i]],useRaster=T,xaxt="n",yaxt="n",ann=FALSE)
    contour(x=sigma_kdens[[i]]$x,y=sigma_kdens[[i]]$y,z=sigma_kdens[[i]]$z,add=T,lty=2)
  }
}

saveGIF(makeplot(),interval=0.2,width=640,height=640)
sigma_col <- colorRampPalette(colors=c("#132B43","#56B1F7"))(50)





makeplot <- function(){
  for(i in 1:length(sigma_kdens)){
    par(bg="#132B43")
    plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(min_x,max_x),ylim=c(min_y,max_y))
    image(sigma_kdens[[i]],useRaster=T,xaxt="n",yaxt="n",ann=FALSE,
          col=sigma_col,xlim=c(min_x,max_x),ylim=c(min_y,max_y),add=TRUE)
    contour(x=sigma_kdens[[i]]$x,y=sigma_kdens[[i]]$y,z=sigma_kdens[[i]]$z,
            add=T,lty=2,col="#9cbfe3",xlim=c(min_x,max_x),ylim=c(min_y,max_y))
  }
}
saveGIF(makeplot(),interval=0.01,ani.width=1024,ani.height=1024)

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

goldpr1 <- adaptMCMC(target=goldpr,init_theta=c(1,1),covmat=diag(c(1,1)),n_iterations=1e3,
             adapt_size_start=10,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=20,
             verbose=FALSE,info=1e2)

par(mfrow=c(1,2))

x1 <- seq(-2, 2, length=100)
x2 <- seq(-2, 2, length=100)
d.goldpr <- matrix(apply(expand.grid(x1, x2), 1, goldpr), nrow=100)
image(x1, x2, d.goldpr, col=cm.colors(100))
contour(x1, x2, d.goldpr, add=TRUE, col=gray(0.6))
lines(goldpr1$theta_trace, type='l')

matplot(goldpr1$theta_trace,type="l")

par(mfrow=c(1,1))
