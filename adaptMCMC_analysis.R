##############################################
######Analysis of Adaptive MCMC routines######
######Sean Wu 8/5/2016########################
##############################################

Rcpp::sourceCpp('/Users/slwu89/Dropbox/GitHub/MCMC/adaptMCMC_source.cpp')

#load packages
library(foreach)
library(doSNOW)
library(parallel)
library(animation)


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


################################################
####Graphics for Rosenbrock (banana) Function###
################################################


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

###Base plotting with Animation package###

max_x <- max(sapply(sigma_kdens,function(x){max(x$x)}))
min_x <- min(sapply(sigma_kdens,function(x){min(x$x)}))
max_y <- max(sapply(sigma_kdens,function(x){max(x$y)}))
min_y <- min(sapply(sigma_kdens,function(x){min(x$y)}))

x1 <- seq(-30,30,length=1e3)
x2 <- seq(-30,30,length=1e3)
x1x2 <- expand.grid(x1,x2)
d.banana <- matrix(apply(x1x2, 1, p.log), nrow=1e3)

sigma_col <- colorRampPalette(colors=c("#132B43","#56B1F7"))(100)

makeplot <- function(){
  for(i in 1:length(sigma_kdens)){
    par(bg="#132B43")
    plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(min_x,max_x),ylim=c(min_y,max_y))
    image(sigma_kdens[[i]],useRaster=T,xaxt="n",yaxt="n",ann=FALSE,
          col=sigma_col,xlim=c(min_x,max_x),ylim=c(min_y,max_y),add=TRUE)
    contour(x=sigma_kdens[[i]]$x,y=sigma_kdens[[i]]$y,z=sigma_kdens[[i]]$z,
            add=T,lty=2,col="#9cbfe3",xlim=c(min_x,max_x),ylim=c(min_y,max_y))
    contour(x1, x2, exp(d.banana), add=TRUE, col="#9cbfe3", drawlabels=FALSE)
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

set.seed(123)
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


##############################################
####Graphics for Goldenstein-Price Function###
##############################################


#only record data when the sampler changed position
goldpr_index <- 1
for(i in 2:nrow(goldpr1$theta_trace)){
  if(any(goldpr1$theta_trace[i,] != goldpr1$theta_trace[i-1,])){
    goldpr_index <- append(goldpr_index,i)
  } else {
    next()
  }
}

#estimate kernel denstiy of empirical sigma
cl <- makeCluster(spec=detectCores()-2)
registerDoSNOW(cl)
sigma_kdens <- foreach(i=goldpr_index,.packages=c("MASS"),.verbose=TRUE) %dopar% {
  if(sum(goldpr1$sigma_empirical[,,i])==0){
    return(NULL)
  }
  reps <- mvrnorm(1e4,mu=goldpr1$theta_trace[i,],Sigma=goldpr1$sigma_empirical[,,i])
  dd <- kde2d(reps[,1],reps[,2],n=200)
  return(dd)
}

stopCluster(cl)
rm(cl)

sigma_kdens <- sigma_kdens[!sapply(sigma_kdens,is.null)]

###preparing data for plotting function contours and sampling distribution###

# max_x <- max(sapply(sigma_kdens,function(x){max(x$x)}))
# min_x <- min(sapply(sigma_kdens,function(x){min(x$x)}))
# max_y <- max(sapply(sigma_kdens,function(x){max(x$y)}))
# min_y <- min(sapply(sigma_kdens,function(x){min(x$y)}))

x1 <- seq(-5.5,5.5,length=500)
x2 <- seq(-5.5,5.5,length=500)
x1x2 <- expand.grid(x1,x2)
d.goldpr <- matrix(apply(x1x2, 1, goldpr), nrow=500)

sigma_col <- colorRampPalette(colors=c("#132B43","#56B1F7"))(100) #colors for sampling distribution
goldpr_col <- heat.colors(50,alpha=0.8) #colors for function contours

###preparing data for plotting trace of chain###

block <- 55 #how many iterations of chain to print
trace_dat <- goldpr1$theta_trace[goldpr_index,]
trace_col <- rev(colorRampPalette(colors=c("#132B43","#56B1F7"),alpha=0.8)(block)) #colors for trace of chain

makeplot <- function(){
  for(i in 1:length(sigma_kdens)){  
    
    #set up background color and remove margins
    par(bg="#132B43",mar=c(0,0,0,0),mgp=c(0,0,0))
    
    #set up empty plot
    plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(-5,5),ylim=c(-5,5))
    
    #plot the sampling distribution raster
    image(sigma_kdens[[i]],useRaster=T,xaxt="n",yaxt="n",ann=FALSE,
          col=sigma_col,xlim=c(-5,5),ylim=c(-5,5),add=TRUE)
    
    #plot the function contours
    contour(x1, x2, d.goldpr, add=TRUE, col=goldpr_col, drawlabels=FALSE,
            nlevels = 50,lwd = 1.25) 
    
    #plot the sampling distribution contours
    contour(x=sigma_kdens[[i]]$x,y=sigma_kdens[[i]]$y,z=sigma_kdens[[i]]$z,
            add=T,lty=2,col="#9cbfe3",xlim=c(-5,5),ylim=c(-5,5))
    
    #plot trace of chain
    if(i > 1){ #can't draw a line of 1 point
      
      if(i <= block){ #beginning block
        current_i <- trace_dat[i:1,] #pull out trace
        for(j in (nrow(current_i)-1):1){ #draw the current segment
          lines(current_i[(j+1):j,],col=trace_col[j])
        }
        
      } else { #moving-window blocks
        
        current_i <- trace_dat[i:(i-(block-1)),] #pull out trace
        for(j in (nrow(current_i)-1):1){ #draw the current segment
          lines(current_i[(j+1):j,],col=trace_col[j])
        }
        
      } 
    } #end plotting trace of chain
    
  } #end for loop
}

saveGIF(makeplot(),interval=0.01,ani.width=1024,ani.height=1024,
        movie.name="goldpr_mcmc.gif")

