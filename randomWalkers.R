#######################################################
######Visualization of C++ Random Walk Algorithms######
######Sean Wu 8/17/2016################################
#######################################################

#souce randomWalkers.cpp

##Random Walk Animation###
library(animation)

random_walker2d <- random_walk2d(5e4)
plot(random_walker2d,type="l",col="blue")

#make_colors will return a smooth interpolation of heat colors to black based on length of current block
make_colors_heat <- function(length){
  colors <- rev(heat.colors(length/2))
  colors <- c(colors,colorRampPalette(c(colors[length(colors)],"black"))(length/2))
  return(colors)
}

make_colors_topo <- function(length){
  colors <- rev(topo.colors(length/2))
  colors <- c(colors,colorRampPalette(c(colors[length(colors)],"black"))(length/2))
  return(colors)
}

make_colors_terrain <- function(length){
  colors <- rev(terrain.colors(length/2))
  colors <- c(colors,colorRampPalette(c(colors[length(colors)],"black"))(length/2))
  return(colors)
}

block <- 200 #how many iterations to plot at once
color_block <- make_colors_heat(block)
end <- 5000

max_x <- max(random_walker2d[1:end,1])
min_x <- min(random_walker2d[1:end,1])
max_y <- max(random_walker2d[1:end,2])
min_y <- min(random_walker2d[1:end,2])

#setup_plot opens a blank plotting surface with the correct boundaries
setup_plot <- function(){
  plot(c(0,0),type="n",axes=FALSE,frame.plot=TRUE,ann=FALSE,xlim=c(min_x,max_x),ylim=c(min_y,max_y))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="black")
}

#generate animation
saveGIF(
  for(i in 2:(end-1)){
    
    if(i <= block){ #beginning block
      
      current_i <- random_walker2d[i:1,]
      setup_plot()
      
      for(j in (nrow(current_i)-1):1){ #draw the current segment
        lines(current_i[(j+1):j,],col=color_block[j])
      }
      
    } else { #moving-window blocks
      
      current_i <- random_walker2d[i:(i-(block-1)),]
      setup_plot()
      
      for(j in (nrow(current_i)-1):1){ #draw the current segment
        lines(current_i[(j+1):j,],col=color_block[j])
      }
    }
    
  }
  ,movie.name="rw2d.gif",ani.height=640,ani.width=640,interval=0.1)

###Animate Multiple Random Walks###
rw1 <- random_walk2d(1e4)
rw2 <- random_walk2d(1e4)
rw3 <- random_walk2d(1e4)

max_x <- max(c(max(rw1[,1]),max(rw2[,1]),max(rw3[,1])))
min_x <- min(c(min(rw1[,1]),min(rw2[,1]),min(rw3[,1])))
max_y <- max(c(max(rw1[,2]),max(rw2[,2]),max(rw3[,2])))
min_y <- min(c(min(rw1[,2]),min(rw2[,2]),min(rw3[,2])))

setup_plot <- function(){
  plot(c(0,0),type="n",axes=FALSE,frame.plot=TRUE,ann=FALSE,xlim=c(min_x,max_x),ylim=c(min_y,max_y))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="black")
}

block <- 750
color_block_heat <- make_colors_heat(block)
color_block_topo <- make_colors_topo(block)
color_block_terrain <- make_colors_terrain(block)
end <- 1e4

#GIF version
saveGIF(
  for(i in 2:(end-1)){
    
    if(i <= block){ #beginning block
      
      current_i_1 <- rw1[i:1,]
      current_i_2 <- rw2[i:1,]
      current_i_3 <- rw3[i:1,]
      setup_plot()
      
      for(j in (nrow(current_i_1)-1):1){ #draw the current segment
        lines(current_i_1[(j+1):j,],col=color_block_heat[j])
        lines(current_i_2[(j+1):j,],col=color_block_topo[j])
        lines(current_i_3[(j+1):j,],col=color_block_terrain[j])
      }
      
    } else { #moving-window blocks
      
      current_i_1 <- rw1[i:(i-(block-1)),]
      current_i_2 <- rw2[i:(i-(block-1)),]
      current_i_2 <- rw3[i:(i-(block-1)),]
      setup_plot()
      
      for(j in (nrow(current_i_1)-1):1){ #draw the current segment
        lines(current_i_1[(j+1):j,],col=color_block_heat[j])
        lines(current_i_2[(j+1):j,],col=color_block_topo[j])
        lines(current_i_3[(j+1):j,],col=color_block_terrain[j])
      }
    }
  }
  ,movie.name="dueling_walkers.gif",ani.height=640,ani.width=640,interval=0.01)

#HTML version
saveHTML({
  for(i in 2:(end-1)){
    
    if(i <= block){ #beginning block
      
      current_i_1 <- rw1[i:1,]
      current_i_2 <- rw2[i:1,]
      current_i_3 <- rw3[i:1,]
      setup_plot()
      
      for(j in (nrow(current_i_1)-1):1){ #draw the current segment
        lines(current_i_1[(j+1):j,],col=color_block_heat[j])
        lines(current_i_2[(j+1):j,],col=color_block_topo[j])
        lines(current_i_3[(j+1):j,],col=color_block_terrain[j])
      }
      
    } else { #moving-window blocks
      
      current_i_1 <- rw1[i:(i-(block-1)),]
      current_i_2 <- rw2[i:(i-(block-1)),]
      current_i_2 <- rw3[i:(i-(block-1)),]
      setup_plot()
      
      for(j in (nrow(current_i_1)-1):1){ #draw the current segment
        lines(current_i_1[(j+1):j,],col=color_block_heat[j])
        lines(current_i_2[(j+1):j,],col=color_block_topo[j])
        lines(current_i_3[(j+1):j,],col=color_block_terrain[j])
      }
    }
  }
},img.name="duel_rw",imgdir="duel_rw_dir",htmlfile="duel_rw.html",
autobrowse=FALSE,title="3 Random Walk Algorithms")

###3D Random Walker###
library(scatterplot3d)
random_walker3d <- random_walk3d(5e4)
scatterplot3d(random_walker3d,type="l",color="blue")


###Random walk on graph###
markov_graph <- matrix(data=0,nrow=100,ncol=100)
for(i in 1:nrow(markov_graph)){
  transitions <- runif(ncol(markov_graph))
  markov_graph[i,] <- transitions/sum(transitions)
}
graph_walker <- random_walkGraph(graph=markov_graph,init_pos=1,n=500)
graph_walker <- graph_walker + 1