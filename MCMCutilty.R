########################################################
######Visualization and Utility Functions for MCMC######
######Sean Wu 4/11/2016#################################
########################################################

###load libraries###
####################

library(reshape2)
library(plyr)
library(coda)
library(ggplot2)
library(parallel)

###visualization functions###
#############################

#convert trace from parallel chains into melted data frame
melt_trace <- function(chain){
  
  #get mcmc trace into data frame
  chain_df <- ldply(.data=chain,.fun=function(x){
    x$trace
  },.progress="text")
  
  #add the chain index and iteration
  chain_df$chain <- unlist(lapply(1:length(chain),function(x){
    rep(x,nrow(chain[[1]]$trace))
  }))
  chain_df$iteration <- rep(1:nrow(chain[[1]]$trace),length(chain))
  
  #return the long format data
  chain_df$chain <- as.factor(chain_df$chain)
  return(chain_df)
}

#view posterior density across multiple chains
univar_posterior <- function(data,par,alpha=0.5){
  
  #create vector of indices for what we need
  index <- as.character(c(par,"iteration","chain"))
  
  #plot the posterior density across all chains
  ggplot(data=data[,index]) +
    geom_density(aes_string(x=par,fill="chain",colour="chain"),alpha=alpha) +
    guides(fill=FALSE,colour=FALSE) +
    labs(x=par,y="Density") +
    theme_bw()
}

#view trace of a parameter across multiple chains
univar_trace <- function(data,par,alpha=0.75){
  
  #create vector of indices for what we need
  index <- as.character(c(par,"iteration","chain"))
  
  #plot the trace across all chains
  ggplot(data=data[,index]) +
    geom_line(aes_string(x="iteration",y=par,colour="chain"),alpha=alpha) +
    guides(colour=FALSE) +
    labs(x="Iteration",y=par) +
    theme_bw()
}

#view covariance matrix of all parameters across multiple chains
multivar_covariance <- function(data,method,low="blue",high="red"){
  
  #calculate covariance matrix 
  cov_mat <- cor(as.matrix(data[,-which(names(data) %in% c("chain","iteration"))]),method=method)
  cov_mat[lower.tri(cov_mat)] <- NA
  cov_mat <- melt(cov_mat[,nrow(cov_mat):1])
  
  ggplot(data=cov_mat) +
    geom_tile(aes(x=Var1,y=Var2,fill=value)) +
    scale_fill_continuous(low=low,high=high,na.value="transparent",guide=guide_colorbar(title=paste(method,"\ncorrelation"),title.position="bottom",barwidth=10)) +
    geom_text(aes(x=Var1,y=Var2,label=round(value,digits=4)),color="white",size=5) +
    theme(legend.title=element_text(size=12),legend.position=c(.70,.75),legend.direction="horizontal",panel.background=element_blank(),panel.border=element_rect(fill=NA,size=1.2,colour="black"),axis.title=element_blank(),axis.text=element_text(size=12))
}

#view bivariate scatterplot of samples from posterior density across multiple chains
scatter_dens <- function(par1,par2,low="blue",high="red",data,...){
  
  #extract ranges of highest density
  init_plot <- ggplot(data=data,aes_string(x=par1,y=par2)) +
    stat_density_2d(aes(fill=..level..,alpha=..level..),geom="polygon",colour="black",...)
  range <- ggplot_build(init_plot)$panel$ranges[[1]][c("x.range","y.range")]
  
  #create plot
  init_plot + scale_fill_continuous(low=low,high=high) +
    guides(alpha=FALSE) +
    coord_cartesian(xlim=c(range[[1]]),ylim=c(range[[2]])) +
    theme_bw() +
    theme(legend.position=c(0,1),legend.justification=c(0,1),legend.background=element_rect(fill="transparent",colour=NA)) +
    labs(fill="Density")
}


###Utility and Helper Functions###
##################################


#parallelsugar mclapply with different random number streams for each socket
mclapply_RNG <- function(
  X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
  mc.silent = FALSE, mc.cores = NULL,
  mc.cleanup = TRUE, mc.allow.recursive = TRUE
) {
  ## Create a cluster
  if (is.null(mc.cores)) {
    mc.cores <- min(length(X), detectCores())
  }
  cl <- parallel::makeCluster( mc.cores )
  
  tryCatch( {
    ## Find out the names of the loaded packages
    loaded.package.names <- c(
      ## Base packages
      sessionInfo()$basePkgs,
      ## Additional packages
      names( sessionInfo()$otherPkgs ))
    
    ### Ship it to the clusters
    parallel::clusterExport(cl,
                            'loaded.package.names',
                            envir=environment())
    
    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parallel::parLapply( cl, 1:length(cl), function(xx){
      lapply(loaded.package.names, function(yy) {
        require(yy , character.only=TRUE)})
    })
    
    parallelsugar:::clusterExport_function(cl, FUN)
    
    #create multiple streams of pseudo-random numbers
    parallel::clusterSetRNGStream(cl=cl,iseed=NULL)
    
    ## Run the lapply in parallel, with a special case for the ... arguments
    if( length( list(...) ) == 0 ) {
      return(parallel::parLapply( cl = cl, X=X, fun=FUN) )
    } else {
      return(parallel::parLapply( cl = cl, X=X, fun=FUN, ...) )
    }
  }, finally = {
    ## Stop the cluster
    parallel::stopCluster(cl)
  })
}
