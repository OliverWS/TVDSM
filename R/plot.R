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


p.format <- function(p){
  if(p < 0.001){
    label = paste0("p ",format.pval(pv=p,digits=3, eps=0.001))
  }
  else {
    label = paste0("p = ",format.pval(pv=p,digits = 3,eps=0.001))
  }
  return(label)
  
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
plotDyad <- function(d, title="",ylabel="EDA",legend_label="Participant", grayscale=F,...) {
  data <- melt(d,id.vars=c("Timestamp"),variable_name = legend_label)
  g <- ggplot(data,aes_string(x="Timestamp",y="value",col=legend_label)) + geom_line(size=1) + ylab(ylabel) + xlab("Time")
  g <- g + ggtitle(title)
  if(grayscale){
    g <- g + scale_color_grey()
  }
  print(g)
  return(g)
}


plotlyDyad <- function(d, title="",ylabel="EDA",legend_label="Participant", grayscale=F,...) {
  data <- melt(d,id.vars=c("Timestamp"),variable_name = legend_label)
  p <- plot_ly(data, x = "Timestamp", y="value",color=legend_label)
  add_lines(p)
  print(p)
  return(p)
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
    
    label = paste0("Shapiro-Wilk Test of Normality\nW = ",format(w$statistic,digits=3), ", ",p.label)
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
    
    y.range <- ggplot_build(plt)$layout$panel_ranges[[1]]$y.range
    x.range <- ggplot_build(plt)$layout$panel_ranges[[1]]$x.range
    
    text.y <- y.range[1] + (y.range[2]-y.range[1])*0.1
    text.x <- x.range[1] + (x.range[2] - x.range[1])*0.5
    plt <- plt + annotate("label",x=text.x,y= text.y,label=try(norm.test(x)),size=5,hjust="center")
  }
  
  
  return(plt)
}

ggbracket <- function(plt,pad = 0.1,offset=0,sig.cutoff=0.05,ci.func=mean_cl_normal){
  f.test <- function(x,y){
    mdl <- lm(y ~ x)
    s <- summary(aov(mdl))
    f = round(s[[1]][[4]][1],digits = 3)
    p.val = s[[1]][[5]][1]
    if(p.val < 0.001){
      p.label = paste0("p ",format.pval(pv=p.val,digits=3, eps=0.001))
    }
    else {
      p.label = paste0("p = ",format.pval(pv=p.val,digits = 3,eps=0.001))
    }
    
    df = paste(s[[1]][[1]],collapse = ",")
    
    str <- paste0("F(",df,") = ",f,", ",p.label)
    return(str)
  }
  #Scale according to plot limits
  y.range <- ggplot_build(plt)$layout$panel_ranges[[1]]$y.range
  
  pad <- 0.04*(y.range[2]-y.range[1])
  
  x <- plt$data[[as.character(plt$mapping$x)]]
  y <- plt$data[[as.character(plt$mapping$y)]]
  groups <- levels(as.factor(x))
  
  x <- simplify2array(x)
  y <- simplify2array(y)
  
  n=0
  g=plt
  
  comparisons <- combn(groups,m=2,simplify = F)
  
  g <- g + annotate("text",x= groups[2], y= max(y.range)+pad*10*length(comparisons),label=f.test(x,y),hjust="center")
  
  for(p in rev(comparisons)){
    group.1 <- p[1]
    group.1.idx <- which(groups == group.1)
    group.2 <- p[2]
    group.2.idx <- which(groups == group.2)
    y.g1 <- subset(y,x==group.1)
    y.g2 <- subset(y,x==group.2)
    
    t <- t.test(y.g1,y.g2)
    
    if(t$p.value < sig.cutoff){
      if(t$p.value < 0.001){
        label = paste0("p ",format.pval(pv=t$p.value,digits=3, eps=0.001))
      }
      else {
        label = paste0("p = ",format.pval(pv=t$p.value,digits = 3,eps=0.001))
      }
      y.1 <- ci.func(y.g1,conf.int=(1-sig.cutoff))$ymax + pad
      y.2 <-  ci.func(y.g2,conf.int=(1-sig.cutoff))$ymax + pad
      y.max <- max(tapply(y,x,FUN = function(a){return(ci.func(a)$ymax)})) + pad*(n+1)*4
      df <- data.frame(x1=c(group.1,group.1,group.2),x2=c(group.1,group.2,group.2),y1=c(y.max-pad,y.max,y.max),y2=c(y.max,y.max,y.max-pad))
      g <- (g + geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=df,col="black"))
      x.label <- (max(c(group.1.idx,group.2.idx)) - min(c(group.1.idx,group.2.idx)))/2.0 + min(c(group.1.idx,group.2.idx))
      g <- g + annotate("text", x = x.label, y = y.max+pad*2, label = label)
      n = n + 1
    }
    
  }
  return(g)
}

# 
# ggbracket <- function(plt,pad = 0.1,offset=0,sig.cutoff=0.05,ci.func=mean_cl_normal){
#   f.test <- function(x,y){
#     mdl <- lm(y ~ x)
#     s <- summary(aov(mdl))
#     f = round(s[[1]][[4]][1],digits = 3)
#     p.val = s[[1]][[5]][1]
#     if(p.val < 0.001){
#       p.label = paste0("p ",format.pval(pv=p.val,digits=3, eps=0.001))
#     }
#     else {
#       p.label = paste0("p = ",format.pval(pv=p.val,digits = 3,eps=0.001))
#     }
#     
#     df = paste(s[[1]][[1]],collapse = ",")
#     
#     str <- paste0("F(",df,") = ",f,", ",p.label)
#     return(str)
#   }
#   #Scale according to plot limits
#   y.range <- ggplot_build(plt)$panel$ranges[[1]]$y.range
#   pad <- 0.05*(y.range[2]-y.range[1])
#   
#   x <- plt$data[[as.character(plt$mapping$x)]]
#   y <- plt$data[[as.character(plt$mapping$y)]]
#   groups <- levels(as.factor(x))
#   
#   x <- simplify2array(x)
#   y <- simplify2array(y)
#   
#   n=0
#   g=plt
#   
#   comparisons <- combn(groups,m=2,simplify = F)
#   
#   g <- g + annotate("text",x= groups[2], y= max(y.range)+pad*10*length(comparisons),label=f.test(x,y),hjust="center")
#   
#   for(p in rev(comparisons)){
#     group.1 <- p[1]
#     group.1.idx <- which(groups == group.1)
#     group.2 <- p[2]
#     group.2.idx <- which(groups == group.2)
#     y.g1 <- subset(y,x==group.1)
#     y.g2 <- subset(y,x==group.2)
#     
#     t <- t.test(y.g1,y.g2)
#     
#     if(t$p.value < sig.cutoff){
#       if(t$p.value < 0.001){
#         label = paste0("p ",format.pval(pv=t$p.value,digits=3, eps=0.001))
#       }
#       else {
#         label = paste0("p = ",format.pval(pv=t$p.value,digits = 3,eps=0.001))
#       }
#       y.1 <- ci.func(y.g1,conf.int=(1-sig.cutoff))$ymax + pad
#       y.2 <-  ci.func(y.g2,conf.int=(1-sig.cutoff))$ymax + pad
#       y.max <- max(tapply(y,x,FUN = function(a){return(ci.func(a)$ymax)})) + pad*(n+1)*5
#       df <- data.frame(x1=c(group.1,group.1,group.2),x2=c(group.1,group.2,group.2),y1=c(y.max-pad,y.max,y.max),y2=c(y.max,y.max,y.max-pad))
#       g <- (g + geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=df,col="black"))
#       x.label <- (max(c(group.1.idx,group.2.idx)) - min(c(group.1.idx,group.2.idx)))/2.0 + min(c(group.1.idx,group.2.idx))
#       g <- g + annotate("text", x = x.label, y = y.max+pad*2, label = label)
#       n = n + 1
#     }
#     
#   }
#   return(g)
# }

plotMetric <- function(data,metric,group,sig.cutoff=0.05) {
  n_fun <- function(x){
    return(data.frame(y = mean(x, na.rm = T) + se(x), label = paste0("n = ",length(x))))
  }
  
  p <- ggplot(data=data, aes_string(group, metric,col=group)) + 
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar",width=0.25,fun.args = list(conf.int=(1-sig.cutoff))) 
  #stat_summary(fun.data = n_fun, geom = "text",vjust="bottom",show.legend = F) + 
  # ggbracket(data[group],data[metric],"MV","TD") + 
  p <- ggbracket(p,sig.cutoff = sig.cutoff,ci.func = mean_cl_normal) 
  # geom_errorbarh(aes(xmin=g1,xmax=g2,y =  as.numeric(g1:g2)*0.25,x=(as.numeric(g2)-as.numeric(g1))/2,col=g1:g2),data=data.frame(g1=c("MV","MV","TD"),g2=c("V","TD","V")),inherit.aes = F,height=0.1,show.legend = F)
  return(p)
  
}

ggrel.x <- function(plt=last_plot(),x=0.5){
  x.range <- ggplot_build(plt)$layout$panel_params[[1]]$x.range
  coord.x = min(x.range) + (max(x.range) - min(x.range))*x
  return(coord.x)
  
}

ggrel.y <- function(plt=last_plot(),y=0.5){
  y.range <- ggplot_build(plt)$layout$panel_params[[1]]$y.range
  coord.y = min(y.range) + (max(y.range) - min(y.range))*y
  return(coord.y)
  
}

annotate.relative <- function(plt=last_plot(),geom="text",x=0.5,y=0.1,...){
  x.range <- ggplot_build(plt)$layout$panel_params[[1]]$x.range
  y.range <- ggplot_build(plt)$layout$panel_params[[1]]$y.range
  coord.x = min(x.range) + (max(x.range) - min(x.range))*x
  coord.y = min(y.range) + (max(y.range) - min(y.range))*y
  
  return(annotate(geom=geom,x = coord.x,y=coord.y,...))
}

scatter_plot <- function(x,y,data){
  cor_fun <- function(x,y,method="pearson"){
    r = NULL
    try((r = cor.test(x,y,method = method)))
    if (is.null(r)) {
      return("r = 0, p = n.s.")
    }
    else {
      return(paste0("r = ",round(r$estimate,2),", ",p.format(r$p.value)))
    }
  }
  
  p<- ggplot(data, aes_string(x=x, y=y)) +
    geom_point() +    # Use hollow circles
    # geom_text(aes(label=ID),nudge_y=0.005) +
    geom_smooth(method="lm",se=T)    # Don't add shaded confidence region
  cor.all <- cor_fun(data[[x]], data[[y]])
  p <- p + annotate.relative(geom="label", plt = p,y=0.95,x=0.9,label=cor.all,hjust="right")
  p
}
