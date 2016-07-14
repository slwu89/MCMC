######Visualization code for mcmc trace######
######Sean Wu 4/11/2016######################
#############################################

###load libraries###
####################

library(reshape2)
library(plyr)
library(coda)
library(ggplot2)

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

#plot results of particle filter (sequential Monte Carlo)
plot_smc <- function(smc,data=NULL,type=NULL,colour="steelblue"){
  
  if(!is.character(type)) {
    stop("Please enter a character of the compartment you want to plot!")
  }

  #prepare data
  traj <- smc$trajectory
  names(traj) <- 1:length(traj)
  
  #extract compartment of interest
  traj_com <- ldply(traj,function(x){
    type_traj <- ddply(x,"time",function(x) x[type])
    names(type_traj)[ncol(type_traj)] <- "V1"
    type_traj <- join(x,type_traj,by="time")
    return(type_traj)
  },.id="replicate",.progress="text")
  
  #compute summary statistics
  traj_CI <- ddply(traj_com,c("time"),function(x){
    tmp_CI <- as.data.frame(t(quantile(x$V1,prob=c(0.025, 0.25, 0.5, 0.75, 0.975))))
    names(tmp_CI) <- c("low_95", "low_50", "median", "up_50", "up_95")
    tmp_CI$mean <- mean(x$V1)
    return(tmp_CI)
  })

  #prepare summary statistics for plotting
  traj_CI_line <- melt(traj_CI[c("time","mean","median")],id.vars="time")
  traj_CI_area <- melt(traj_CI[c("time","low_95","low_50","up_50","up_95")],id.vars="time")
  traj_CI_area$type <- sapply(traj_CI_area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj_CI_area$CI <- sapply(traj_CI_area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj_CI_area$variable <- NULL
  traj_CI_area <- dcast(traj_CI_area, paste0("time","+CI~type"))
  
  #plot summary statistics
  p <- ggplot(traj_CI_area)
  p <- p + geom_ribbon(data = traj_CI_area, aes_string(x = "time", ymin = "low", ymax = "up", alpha = "CI"), fill = colour)
  p <- p + geom_line(data = traj_CI_line, aes_string(x = "time", y = "value", linetype = "variable"), colour = colour)
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.25, "50" = 0.45), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  #add data (if present)
  if(!is.null(data)){
    print("Adding data points")
    p <- p + geom_point(data=data,aes(x=time,y=I),colour="black")
  }
  
  #print to screen
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  print(p)
}