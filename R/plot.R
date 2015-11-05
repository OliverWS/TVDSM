require(ggplot2)
library(grid)
library(signal)
library(gridExtra)

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
      acc <- ggplot(data=data) + geom_line(aes(x = Timestamp, y=X),colour="red") + geom_line(aes(x = Timestamp, y=Y),colour="green") + geom_line(aes(x = Timestamp, y=Z),color="blue") + xlab("Time") + ylab("Acc") +  ggtitle("Accelerometer")
      
    }
    else{
      data$Movement <- sqrt(data$X*data$X+data$Y*data$Y+data$Z*data$Z)
      acc <- ggplot(data=data) + geom_line(aes(x = Timestamp, y=Movement),colour="red") + xlab("Time") + ylab("Movment") +  ggtitle("Accelerometer")
      
    }
    return(multiplot(eda,acc))
    
  }
  else {
    return(eda)
  }
  
  
  
}


plot.dyad <- function(p1,p2,title="Dyad") {
  colors <- c("#2CC4FF","#DD6DAE")
  p1.name <- deparse(substitute(p1))
  p2.name <- deparse(substitute(p2))
  
  
  data <- as.dyad(p1,p2)
  cols <- colnames(data)
  colnames(data) <- c("Timestamp","P1.EDA","P2.EDA")
  ggplot(data=data, aes(x=Timestamp, y=P1.EDA)) + geom_line(aes(x = Timestamp, y=P1.EDA),colour="#2CC4FF")+ geom_line(aes(x = Timestamp, y=P2.EDA),colour="#DD6DAE")+ xlab("Time") + ggtitle(paste(p1.name,"&",p2.name)) + ylab("EDA (uS)")
}


plotEDAByCondition <- function(eda, codes, title="") {
  
  if(is.character(eda)) {
    eda <- read.eda(eda)
  }

  g1 <- ggplot(data=eda, aes(x=Timestamp, y=EDA)) + geom_line(aes(x = Timestamp, y=EDA),colour="#2CC4FF")+  xlab("Time") + ggtitle(title) + ylab("EDA (uS)")
  g1 <- g1 + geom_rect(aes(xmin=`Start Time`,xmax=`End Time`,ymin=-1.0,ymax=0, fill=Condition,col=NULL),data=codes)
  return(g1)
  
}