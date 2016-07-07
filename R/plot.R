library(ggplot2)
library(grid)
library(signal)
library(gridExtra)
library(reshape)
library(cowplot)
# Multiple plot function
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#



plot.ci <- function(formula){ 



}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,heights=c(2,1)) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout),heights = heights)))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot.eda <- function(data,title="EDA Data",show_acc=F, show_raw_acc=F) {
    
  eda <- ggplot(data=data, aes(x=Timestamp,y=EDA)) + geom_line(aes(x = Timestamp, y=EDA), colour="#2CC4FF") + xlab("Time") + ylab("EDA") + ggtitle(title) #+ geom_area(aes(x = Timestamp, y=EDA), fill="#B6E1F2")
  if(show_acc) {
    if(show_raw_acc){
      acc <- ggplot(data=data) + geom_line(aes(x = Timestamp, y=X),colour="red") + geom_line(aes(x = Timestamp, y=Y),colour="green") + geom_line(aes(x = Timestamp, y=Z),color="blue") + xlab("Time") + ylab("Acc") 
      
    }
    else{
      data$Movement <- sqrt(data$X*data$X+data$Y*data$Y+data$Z*data$Z)
      acc <- ggplot(data=data) + geom_line(aes(x = Timestamp, y=Movement),colour="red") + xlab("Time") + ylab("Motion")
      
    }
    return(plot_grid(eda,acc,ncol=1,rel_heights = c(3,1)))
    
  }
  else {
    return(eda)
  }
  
  
  
}
plotDyad <- function(d, title="",ylabel="EDA",legend_label="Participant") {
  data <- melt(d,id.vars=c("Timestamp"),variable_name = legend_label)
  g <- ggplot(data,aes_string(x="Timestamp",y="value",col=legend_label)) + geom_line(size=1) + ylab(ylabel) + xlab("Time")
  g <- g + ggtitle(title)
  print(g)
  return(g)
}





plotEDAByCondition <- function(eda, codes, title="") {
  
  eda <- read.eda(eda)
  g1 <- ggplot(data=eda, aes(x=Timestamp, y=EDA)) + geom_line(aes(x = Timestamp, y=EDA),colour="#2CC4FF")+  xlab("Time") + ggtitle(title) + ylab("EDA (uS)")
  g1 <- g1 + geom_rect(data=codes, aes(xmin=`Start Time`,xmax=`End Time`,ymin=-1.0,ymax=0, fill=Condition,col=NULL))
  
  print(g1)
  return(g1)
  
}


gghist <- function(x,groups=NULL) {
  
  norm.test <- function(data) {
    w <- shapiro.test(data)
    if(w$p.value < 0.001){
      p.label = paste0("p ",format.pval(pv=w$p.value,digits=3, eps=0.001))
    }
    else {
      p.label = paste0("p = ",format.pval(pv=w$p.value,digits = 3,eps=0.001))
    }
    
    label = paste0("W = ",format(w$statistic,digits=3), ", ",p.label)
   return( label)
  }
  x.mu <- mean(x,na.rm=T)
  x.sd <- sd(x,na.rm=T)
  if(!is.null(groups)) {
    plt <- ggplot(aes(x=x,fill=group),data = data.frame(x=x,group=groups)) + geom_density(alpha=0.5)

  }
  else {
    plt <- ggplot(aes(x=x),data = data.frame(x=x)) + geom_density(col="darkgray",fill="gray",lwd=1)
    plt <- plt + stat_function(fun = dnorm, colour = "red",args = list(mean=x.mu,sd=x.sd))
    
    y.range <- ggplot_build(plt)$panel$ranges[[1]]$y.range
    x.range <- ggplot_build(plt)$panel$ranges[[1]]$x.range
    
    text.y <- y.range[1] + (y.range[2]-y.range[1])*0.1
    text.x <- x.range[1] + (x.range[2] - x.range[1])*0.5
    plt <- plt + annotate("text",x=text.x,y= text.y,label=norm.test(x),size=5)
  }
  
  
  return(plt)
}
