library(lattice)
library(lme4)
library(psychometric)
library(reshape)
library(reshape2)
library(e1071)
library(lmerTest)
library(ggplot2)
library(plyr)
library(cowplot)
library(lmtest)
library(reshape)
library(systemfit)
library(Hmisc)
library(progress)

#' Run ANOVA on model set
#'
#' @description Runs ANOVA comparisons on delta R-squared values from a list of models using a specified condition.
#' @param mdls A list of model objects.
#' @param key The key to extract the condition variable (default "Condition").
#' @return A fitted lmer model object.
#' @export

runAnova <- function(mdls, key="Condition"){
  condition <- get.key(mdls,key)
  
  b5 <- get.key(mdls,"b5")
  b2<- get.key(mdls,"b2")
  b21 <- get.key(mdls,"b21")
  b45 <- get.key(mdls,"b45")

  x.dr2 <- get.key(mdls,"x.r.squared") - get.key(mdls,"x.base.r.squared")
  y.dr2 <- get.key(mdls,"y.r.squared") - get.key(mdls,"y.base.r.squared")
  time <- get.key(mdls,key="Timestamp")
  dyad <- get.key(mdls,"dyad.name")
  t <- seq_along(x.dr2)
  condition <- as.factor(condition)
  print(table(condition))
  print(summary(lm(x.dr2 ~ condition)),digits=3)
  
  #   mmdl.1 <- lmer(avg.dr2 ~ SJ + Ob + (1 | dyad),REML=FALSE)
  #   print(summary(mmdl.1))
  mmdl.x.empty <- lmer(x.dr2 ~ (1 | dyad),REML=F) 
  mmdl.x <- lmer(x.dr2 ~ condition + (1 | dyad),REML=F) 
  print(summary(mmdl.x))
  print(anova(mmdl.x.empty,mmdl.x))
  
  mmdl.y.empty <- lmer(y.dr2 ~ (1 | dyad),REML=F) 
  mmdl.y <- lmer(y.dr2 ~ condition + (1 | dyad),REML=F) 
  print(summary(mmdl.y))
  print(anova(mmdl.y.empty,mmdl.y))
  
  return(mmdl.x)
}

#' Validate condition times within data range
#'
#' @description Checks each condition's start and end times to confirm they fall within the dyad's timestamp range.
#' @param d A data frame with a Timestamp column.
#' @param codes A data frame with condition timing (with Start.Time, End.Time, and Condition).
#' @param timeformat A string specifying the timestamp format (default "%Y-%m-%d %H:%M:%S").
#' @return Stops with an error if a condition time is out-of-range.
#' @export

validateConditionTimes <- function(d,codes,timeformat ="%Y-%m-%d %H:%M:%S"){
  start.ctime <- as.POSIXct(d$Timestamp[[1]])
  end.ctime <- as.POSIXct(d$Timestamp[[length(d$Timestamp)]])
  
  codes$Start.Time.CTime <- as.POSIXct(codes$Start.Time,format=timeformat)
  codes$End.Time.CTime <- as.POSIXct(codes$End.Time,format=timeformat)
  
  for(n in 1:dim(codes)[[1]]){
    st = codes$Start.Time.CTime[[n]]
    et = codes$End.Time.CTime[[n]]
    if((st < start.ctime) | (et > end.ctime)){
      msg <- paste0("Error: times for condition '",codes$Condition[[n]],"' (",strftime(st,format = timeformat), " to ", strftime(et,format=timeformat),") are outside start/end times for dyad (",strftime(start.ctime,format=timeformat)," to ",strftime(end.ctime,format=timeformat),")")
      stop(msg,call. = F,domain = NA)
    }
  }
}

#' Extract a key from a list of objects
#'
#' @description Retrieves the specified key from each element of a list.
#' @param lst A list of objects.
#' @param key The key to extract.
#' @param as.list Logical; if TRUE, returns a list (default FALSE).
#' @return A vector or list of values corresponding to the key from each element.
#' @export

get.key <- function(lst, key,as.list=F) {
  if(as.list == T){
    return(as.list(lapply(lst,FUN = function(l){return(try(l[[key]]))}),recursive = F ))
  }
  else{
    return(unlist( lapply(lst,FUN = function(l){return(try(l[[key]]))}),recursive = F ))
  }
}

#' Plot lag parameters over time
#'
#' @description Plots the lags and R-squared values over time using ggplot2.
#' @param data A data frame containing timestamps and parameter values.
#' @param xname Label for the x variable.
#' @param yname Label for the y variable.
#' @param title The plot title.
#' @param save Logical, whether to save the plot.
#' @param f Filename to save the plot (default "plot.pdf").
#' @param use.delta.rsquared Logical, whether to use delta R-squared values.
#' @param by.condition Logical, whether to group by condition.
#' @param autoscale Logical, whether to auto-scale the plot.
#' @return A combined ggplot object.
#' @export

plot.lagparams <- function(data,xname="x",yname="y",title=paste(xname,"&",yname),save=F,f="plot.pdf",use.delta.rsquared=T,by.condition=T, autoscale=T) {
  data <- na.omit(data)
  r2.labels <- c(bquote(I[.(xname)%<-% .(yname)]^2),bquote(I[.(yname)%<-%.(xname)]^2))
  data$timestamp <- as.POSIXct(data$timestamp,origin = "1970-01-01")
  params1 <- melt(subset(data,select=c("timestamp","x.lag","y.lag")),id.vars = c("timestamp"))
  params2 <- melt(subset(data,select=c("timestamp","dx.r.squared","dy.r.squared")),id.vars = c("timestamp"))

  
  glegend1 <- scale_colour_discrete(name  ="Variable",
                                    breaks=c("x.lag","y.lag"),
                                    labels=c(bquote(Lag[.(xname)]),bquote(Lag[.(yname)])))
  glegend2 <- scale_colour_discrete(name  ="Variable",
                                    breaks=c("dx.r.squared","dy.r.squared"),
                                    labels=r2.labels)
  gstyle <- guides(colour = guide_legend(
    label.theme = element_text(size=15,angle=0),label.hjust=0,label.vjust=0),title.theme=element_text(size=15,angle=0))
  x_min = as.POSIXct(min(data$timestamp),origin = "1970-01-01")
  x_max = as.POSIXct(max(data$timestamp),origin = "1970-01-01")
  
  
  plt1 <- ggplot(data=params1, aes(x=timestamp,y=value,colour=variable)) +  geom_line(size=3) + geom_point(size=6) + xlab("Time") + ylab("Lag") + ggtitle(title) +glegend1 +gstyle 
  plt2 <- ggplot(data=params2, aes(x=timestamp,y=value,colour=variable)) + geom_line(size=1) + geom_point(size=3) + xlab("Time") + ylab(bquote(R^2)) + ggtitle(title) + glegend2 + gstyle

  plt1 <- plt1 + scale_x_datetime(limits=c(x_min, x_max))
  plt2 <- plt2 + scale_x_datetime(limits=c(x_min, x_max)) 
  combinedPlt <- plot_grid(plt2,plt1,align = "hv",ncol = 1,nrow = 2,rel_heights = c(2,1))

  if(save){
    ggsave(f)
  }
  return(combinedPlt)
  
}

#' Plot TVDSM model results
#'
#' @description Generates a combined plot displaying TVDSM model parameters and R-squared values over time.
#' @param data A list of TVDSM model objects.
#' @param xname Friendly name for participant 1.
#' @param yname Friendly name for participant 2.
#' @param title Plot title.
#' @param save Logical, whether to save the plot.
#' @param f Filename to save the plot.
#' @param use.delta.rsquared Logical, whether to use delta R-squared values.
#' @param by.condition Logical, whether to group by condition.
#' @param autoscale Logical, whether to auto-scale the plot.
#' @param plotParams Logical, whether to include model parameter plots.
#' @param grayscale Logical, whether to use a grayscale color scheme.
#' @param rawData Optional raw data for plotting.
#' @param returnSeparatePlots Logical, if TRUE returns separate plots.
#' @param fontFamily Font family for text.
#' @return A ggplot object or a list of ggplot objects.
#' @export

plot.tvdsm <- function(data,xname="x",yname="y",title=paste(xname,"&",yname),save=F,f="plot.pdf",use.delta.rsquared=T,by.condition=T, autoscale=F,plotParams=T, grayscale=F,rawData=NULL,returnSeparatePlots=F,fontFamily="Helvetica") {
  data <- na.omit(data)
  r2.labels <- c(bquote(R[.(xname)]^2),bquote(R[.(yname)]^2))
  if(use.delta.rsquared){
    x.r2 <- get.key(data, key="dx.r.squared")
    y.r2 <- get.key(data, key="dy.r.squared")
    r2.labels <- c(bquote(I[.(xname)%<-% .(yname)]^2),bquote(I[.(yname)%<-%.(xname)]^2))
  }
  else {
    x.r2 <- get.key(data, key="x.r.squared")
    y.r2 <- get.key(data, key="y.r.squared")
  }
  b2 <- get.key(data, key="b2")
  b5 <- get.key(data, key="b5")
  b1 <- get.key(data, key="b1")
  b4 <- get.key(data, key="b4")
  b45 <- get.key(data, key="b45")
  b21 <- get.key(data, key="b21")
  if(by.condition){
    ES <- get.key(data, key="Condition")
  }
  else {
    ES <- seq_along(b2)
  }
  timestamps <- as.POSIXct(get.key(data,key="Timestamp"),origin = "1970-01-01")
  start <- as.POSIXct(get.key(data,key="Start"),origin = "1970-01-01")
  end <- as.POSIXct(get.key(data,key="End"),origin = "1970-01-01")
  
  params.df <- data.frame(Timestamps=timestamps,b1=b1,b2=b2,b4=b4,b5=b5,b21=b21,b45=b45,x.r2=x.r2,y.r2=y.r2,Condition=ES,start=start,end=end)
  
  len <- length(params.df$Timestamps)
  params1 <- melt(subset(params.df,select=c("Timestamps","b1","b2", "b4","b5","b21","b45")),id.vars = c("Timestamps"))
  params2 <- melt(subset(params.df,select=c("Timestamps","x.r2","y.r2","Condition","start","end")),id.vars = c("Timestamps","start","end","Condition"))
  
  if(all(is.na(params.df$b21)) & all(is.na(params.df$b45))){
    param.labels <- c(bquote(.(xname)[S-R]), bquote(.(xname)[Co-R]),bquote(.(yname)[S-R]),bquote(.(yname)[Co-R]))
    glegend1 <- scale_colour_discrete(name  = "Beta Coefficents",
                                      breaks=c("b1","b2", "b4","b5"),
                                      labels=param.labels)
    glegend1.gray <- scale_shape_discrete(name  = "Beta Coefficents",
                                        breaks=c("b1","b2", "b4","b5"),
                                        labels=param.labels)
    
    
  }
  else {
    param.labels <- c(bquote(.(xname)[S-R]), bquote(.(xname)[Co-R]), bquote(.(xname)[I]), bquote(.(yname)[S-R]),bquote(.(yname)[Co-R]),bquote(.(yname)[I]))
    
  glegend1 <- scale_colour_discrete(name  = "Coefficients",
                                      breaks=c("b1","b2","b21", "b4","b5","b45"),
                                      labels=param.labels)
  glegend1.gray <-  scale_shape_discrete(name  = "Coefficients",
                                                 breaks=c("b1","b2","b21", "b4","b5","b45"),
                                                 labels=param.labels)
    
  
  }
  glegend2 <- scale_colour_discrete(name  = "Effect Size",
                                    breaks=c("x.r2","y.r2"),
                                    labels=r2.labels)
  glegend2.gray <-scale_shape_discrete(name  = "Effect Size",
                                     breaks=c("x.r2","y.r2"),
                                     labels=r2.labels)
    
  

  if(grayscale){
    gstyle <- guides(colour = guide_legend(override.aes=list(fill=NA,size=5,linetype=0),
                                           label.theme = element_text(size=15,angle=0),label.hjust=0,label.vjust=0),title.theme=element_text(size=15,angle=0))
  }
  else {
    gstyle <- guides(colour = guide_legend(override.aes=list(fill=NA),
                                           label.theme = element_text(size=15,angle=0),label.hjust=0,label.vjust=0),title.theme=element_text(size=15,angle=0))
  }

  
  x_min = min(start)
  x_max = max(end)
  point_size = 3*(50.0/len)
  if(point_size > 3) point_size = 3.0
  x_label <- "Time"
  title.1 <- title
  if(plotParams == T){
    title.2 <- NULL
  }
  else if(plotParams == F){
    title.2 <- title
  }
  else {
    title.2 <- NULL
  }
  
  line_size = 1
  if(grayscale) line_size = 0.75
  plt1 <- ggplot(data=params1, aes(x=Timestamps,y=value,colour=variable)) + geom_line(size=line_size) + xlab(x_label) + ylab(NULL) + ggtitle(title.1) +glegend1 +gstyle 
  
  plt2 <- ggplot(data=params2, aes(x=Timestamps,y=value,colour=variable)) + geom_line(size=line_size) + xlab(x_label) + ylab(bquote(I^2)) + ggtitle(title.2) + glegend2 + gstyle
  
  if(point_size > 1.25) {
    if(grayscale){
      plt1 <- plt1 + geom_point(size=point_size,aes(shape=variable)) 
      plt2 <- plt2 + geom_point(size=point_size,aes(shape=variable)) 
    }
    else {
      plt1 <- plt1 + geom_point(size=point_size) 
      plt2 <- plt2 + geom_point(size=point_size) 
    }
  }
  plt1 <- plt1 + scale_x_datetime(limits=c(x_min, x_max)) + theme(text=element_text(family=fontFamily))
  plt2 <- plt2 + scale_x_datetime(limits=c(x_min, x_max)) + theme(text=element_text(family=fontFamily),axis.title.y = element_text(angle = 0))
  
  if(grayscale){
    plt1 <- plt1 + glegend1.gray
    plt2 <- plt2 + glegend2.gray
  }
  
  if(by.condition){
    CONDITION_HEIGHT <- ggrel.y(plt=plt2, y= 0.2)
    if(is.nan(CONDITION_HEIGHT)){
      CONDITION_HEIGHT = 0.2
    }
    
    y.min <- ggrel.y(plt=plt2, y= 0)
    if(is.nan(y.min)) {
      y.min <- 0
    }
    
    plt2 <- plt2 + geom_rect(aes(fill=Condition,ymin=(y.min - CONDITION_HEIGHT),ymax=y.min,xmin=start,xmax=end),linetype=0,alpha=0.8) 
    if(grayscale){
      plt2 <- plt2 + scale_fill_grey()
    }
    if(autoscale == F){
      plt2 <- plt2 + ylim(ggrel.y(plt=plt2, y= 0),1)
    }
    
  }
  else {
    if(autoscale == F){
      plt2 <- plt2 + ylim(0,1)
    }
  }
  if(plotParams == T){
    combinedPlt <- plot_grid(plt1,plt2,align = "hv",ncol = 1,nrow = 2,rel_heights = c(2,1))
  }
  else if(plotParams == "raw"){
    colnames(rawData) <- c("Timestamp",xname,yname)
    plt.raw <- plotDyad(rawData,title = title,ylabel = "EDA (µS)")
    combinedPlt <- plot_grid(plt.raw,plt2,align = "hv",ncol = 1,nrow = 2,rel_heights = c(2,1))
  }
  else {
    combinedPlt <- plt2
  }

  if(save){
    ggsave(f)
  }
  
  if(returnSeparatePlots) {
    return(list(I2Plt=plt2,paramPlt=plt1))
  }
  
  else {
    return(combinedPlt)
  }
  
}

#' Plot density contour of beta coefficients
#'
#' @description Creates a contour plot from a 2D kernel density estimate of beta coefficients.
#' @param data A data frame containing beta coefficient values.
#' @param group Group label for the plot.
#' @param xname Label for the x-axis (default: "Particpant 1 Beta").
#' @param yname Label for the y-axis (default: "Participant 2 Beta").
#' @return A plot of the density contours.
#' @export

plot.density <- function(data,group="Dyad", xname='Particpant 1 Beta', yname="Participant 2 Beta") {
  x.beta <- get.key(data, key="b2")
  y.beta <- get.key(data, key="b4")

  kp1 <- kde2d(x.beta,y.beta,n=500)
  contour(kp1, col='blue',)
  title(paste("Kernel Density Plots -", group), ylab=yname, xlab=xname)
  
  gg <- ggplot(data=kp1,mapping = aes(x,y,z)) + geom_contour()
  print(gg)
}

#' Fit a dynamic system model
#'
#' @description Fits a dynamic system model to time series data using systemfit and SUR.
#' @param x Numeric vector representing x values.
#' @param y Numeric vector representing y values.
#' @param x_mu Baseline value for x; default is mean(x, na.rm=TRUE).
#' @param y_mu Baseline value for y; default is mean(y, na.rm=TRUE).
#' @param type Model type (default 2).
#' @param step Step size parameter (default 0.25).
#' @param p.value Significance level for the likelihood ratio test.
#' @param verbose Logical, whether to print verbose output.
#' @param lag Lag value to apply.
#' @return A list containing the fitted models, R-squared values, and model coefficients.
#' @export

fitDSModel <- function(x,y,x_mu=mean(x,na.rm=T),y_mu=mean(y,na.rm=T),type=2, step=0.25,p.value=0.01,verbose=F,lag=0){
  if(is.null(lag)){
    x_prime <-o.Lag(x,lag=0)
    y_prime <- o.Lag(y,lag=0)
    x_prime.zLag <- o.Lag(x,lag=0)
    y_prime.zLag <- o.Lag(y,lag=0)
    
  }
  else {
    x_prime <-o.Lag(x,lag=lag)
    y_prime <- o.Lag(y,lag=lag)
    x_prime.zLag <- o.Lag(x,lag=0)
    y_prime.zLag <- o.Lag(y,lag=0)
    
  }

  s1 <- ""
  s2 <- ""
  
  data <- data.frame(x_prime=x_prime,y_prime=y_prime, x=x,y=y,x_prime.zLag=x_prime.zLag,y_prime.zLag=y_prime.zLag)
  data <- na.omit(data)
  base_x_model <- x_prime ~ I(x_mu-x)
  base_y_model <- y_prime ~ I(y_mu-y)
  
  
  base_eq <- list(base_x_model,base_y_model)
  basefit <- systemfit(base_eq,data=data, method = "SUR")
  
  base_x_model <- basefit$eq[[1]]
  base_y_model <- basefit$eq[[2]]
  if (type==5){
    
    base_x_model <- x_prime ~ 0 + I(x_mu-x)
    base_y_model <- y_prime ~ 0 + I(y_mu-y)
    
    
    base_eq <- list(base_x_model,base_y_model)
    basefit <- systemfit(base_eq,data=data, method = "SUR")
    
    base_x_model <- basefit$eq[[1]]
    base_y_model <- basefit$eq[[2]]
    
    x_model <- x_prime ~ 0 + I(x_mu-x) + I(y-x)
    y_model <- y_prime ~ 0 + I(y_mu-y) + I(x-y)
    
    
    eq <- list(x_model,y_model)
    fit <- systemfit(eq,data=data, method = "SUR")
    
    x_model <- fit$eq[[1]]
    y_model <- fit$eq[[2]]
    
    b0 <- NA
    b1 <- o.coef(x_model,1)
    b2 <- o.coef(x_model,2)
    b21 <- NA
    b3 <- NA
    b4 <- o.coef(y_model,1)
    b5 <- o.coef(y_model,2)
    b45 <- NA
  }
  
  else if (type==4){
    
    x_model <- x_prime ~ I(x_mu-x) + y_prime.zLag
    y_model <- y_prime ~ I(y_mu-y) + x_prime.zLag
    
    
    eq <- list(x_model,y_model)
    fit <- systemfit(list(x_model,y_model),data=data, method = "SUR")
    
    x_model <- fit$eq[[1]]
    y_model <- fit$eq[[2]]
    
    b0 <- o.coef(x_model,1)
    b1 <- o.coef(x_model,2)
    b2 <- o.coef(x_model,3)
    b21 <- NA
    b3 <- o.coef(y_model,1)
    b4 <- o.coef(y_model,2)
    b5 <- o.coef(y_model,3)
    b45 <- NA
  }
  else if (type==4.5){
    
    x_model <- x_prime ~ I(x_mu-x) * y_prime
    y_model <- y_prime ~ I(y_mu-y) * x_prime
    
    
    eq <- list(x_model,y_model)
    fit <- systemfit(list(x_model,y_model),data=data, method = "SUR")
    
    x_model <- fit$eq[[1]]
    y_model <- fit$eq[[2]]
    
    b0 <- o.coef(x_model,1)
    b1 <- o.coef(x_model,2)
    b2 <- o.coef(x_model,3)
    b21 <- o.coef(x_model,4)
    b3 <- o.coef(y_model,1)
    b4 <- o.coef(y_model,2)
    b5 <- o.coef(y_model,3)
    b45 <- o.coef(y_model,4)
  }
  
  
  else if (type==3){
    
    x_model <- x_prime ~ I(x_mu-x) + I(y-x)
    y_model <- y_prime ~ I(y_mu-y) + I(x-y)
    
    
    eq <- list(x_model,y_model)
    fit <- systemfit(eq,data=data, method = "SUR")
    
    x_model <- fit$eq[[1]]
    y_model <- fit$eq[[2]]
    
    b0 <- o.coef(x_model,1)
    b1 <- o.coef(x_model,2)
    b2 <- o.coef(x_model,3)
    b21 <- NA
    b3 <- o.coef(y_model,1)
    b4 <- o.coef(y_model,2)
    b5 <- o.coef(y_model,3)
    b45 <- NA
  }
  else if(type==2){
    
    x_model <- x_prime ~ I(x_mu-x) * I(y-x)
    y_model <- y_prime ~ I(y_mu-y) * I(x-y)
    
    
    eq <- list(x_model,y_model)
    fit <- systemfit(list(x_model,y_model),data=data, method = "SUR")
    
    x_model <- fit$eq[[1]]
    y_model <- fit$eq[[2]]
    
    b0 <- o.coef(x_model,1)
    b1 <- o.coef(x_model,2)
    b2 <- o.coef(x_model,3)
    b21 <- o.coef(x_model,4)
    b3 <- o.coef(y_model,1)
    b4 <- o.coef(y_model,2)
    b5 <- o.coef(y_model,3)
    b45 <- o.coef(y_model,4)
  }
  
  else {
    x_model <- x_prime ~  I(y-x) 
    y_model <- y_prime ~  I(x-y)
    
    
    eq <- list(x_model,y_model)
    fit <- systemfit(list(x_model,y_model),data=data, method = "SUR")
    
    x_model <- fit$eq[[1]]
    y_model <- fit$eq[[2]]
    print(summary(x_model))
    print(summary(y_model))
    
    
    b0 <- o.coef(x_model,1)
    b1 <- 0
    b2 <- o.coef(x_model,2)
    b3 <- o.coef(y_model,1)
    b4 <- 0
    b5 <- o.coef(y_model,2)
    b21 <- 0
    b45 <- 0
  }
  
 
  if(verbose){
    
    print(summary(x_model))
    print(summary(y_model))
    
    
  }
  tst <- lrtest(basefit,fit)
  if(tst$`Pr(>Chisq)`[[2]] < p.value) {
    x.base.r.squared=summary(base_x_model)$r.squared
    y.base.r.squared=summary(base_y_model)$r.squared
    x.r.squared = summary(x_model)$r.squared
    y.r.squared = summary(y_model)$r.squared
  }
  else {
    if(verbose){
      print(tst)
    }
    x.base.r.squared=0
    y.base.r.squared=0
    
    x.r.squared = 0
    y.r.squared = 0
  }
  
  
  return(list(x_model=x_model,y_model=y_model,base_x_model=base_x_model,base_y_model=base_y_model,model=fit,basemodel=basefit,x.base.r.squared=x.base.r.squared,y.base.r.squared=y.base.r.squared, x.r.squared=x.r.squared, y.r.squared=y.r.squared,dx.r.squared=(x.r.squared - x.base.r.squared),dy.r.squared=(y.r.squared - y.base.r.squared),x.eq=s1,y.eq=s2,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,b5=b5,b21=b21,b45=b45))
}

#' Fit dynamic system model for dyad data
#'
#' @description Downsamples dyad data and fits a dynamic system model for both participants.
#' @param dyad A data frame with dyad time series data.
#' @param type Model type (default 4).
#' @param downsample Downsampling factor (default 1).
#' @param lag Lag value (default 0).
#' @param x_mu Baseline for participant 1.
#' @param y_mu Baseline for participant 2.
#' @param verbose Logical, whether to output verbose messages.
#' @return A list containing model outputs and timing information.
#' @export

fitDSModelForDyad <- function(dyad,type=4,downsample=1,lag=0,x_mu=NULL,y_mu=NULL,verbose=F) {
  FS <- round(getFS(dyad))
  ax <- decimate(dyad[,2], downsample*FS)
  ay <- decimate(dyad[,3],downsample*FS)
  newFS <- round(1.0/downsample)
  mdl <- list(
    x.base.r.squared=0,
    y.base.r.squared=0,
    x.r.squared=0, 
    y.r.squared=0,
    dx.r.squared=NA,
    dy.r.squared=NA,
    x.eq="s1",y.eq="s2",
    b0=NA,b1=NA,b2=NA,b3=NA,b4=NA,b5=NA,b21=NA,b45=NA)
  
tryCatch(mdl <- fitDSModel(ax,ay,p.value = 0.05,type=type,lag=lag*newFS,x_mu = x_mu,y_mu=y_mu,verbose=verbose),
           error=function(e){
             print(e)
           })
  mdl$Timestamp <- min(dyad$Timestamp)
  mdl$Duration <- (max(dyad$Timestamp) - min(dyad$Timestamp))
  mdl$Start <- min(dyad$Timestamp)
  mdl$End <- max(dyad$Timestamp)

  return(mdl)
}

#' Generate URL for model visualization
#'
#' @description Constructs a local file URL for visualizing model results based on coefficients.
#' @param m A model object containing coefficients and a matrix.
#' @return A string representing the URL for the model visualization.
#' @export

urlForModel <- function(m) {
  a1 <- m$b1
  a2 <- m$b2
  b1 <- m$b4
  b2 <- m$b5
  minV <- min(m$mat[,1:2],na.rm = T)
  maxV <- max(m$mat[,1:2], na.rm = T)
  return(paste("file:///Users/OliverWS/git/DCS%20Visualization/index.html?a1=", a1,"&a2=",a2,"&b1=",b1,"&b2=",b2,"&min=",minV,"&max=",maxV,sep=""))
}

#' Analyze model performance over varying lags
#'
#' @description Iterates over a range of lag values to fit models and determine optimal lag based on R-squared.
#' @param f1 File path for participant 1.
#' @param f2 File path for participant 2.
#' @param dyad Optional dyad object (if provided, f1 and f2 are ignored).
#' @param xname Friendly name for participant 1.
#' @param yname Friendly name for participant 2.
#' @param norm Logical, whether to standardize data.
#' @param window_size Window size in seconds.
#' @param window_step Step size in seconds.
#' @param start Start time.
#' @param end End time.
#' @param func Function to use for fitting the model.
#' @param na.rm Logical, whether to remove NA values.
#' @param simulate Logical, whether to simulate dyad data.
#' @param dname Friendly name for the dyad.
#' @param measure Measurement column (e.g., "EDA").
#' @param downsample Downsampling factor.
#' @param minLag Minimum lag value.
#' @param maxLag Maximum lag value.
#' @param noPlots Logical, if TRUE no plots are generated.
#' @param relativeToLag Reference lag index (default -1).
#' @param type Model type.
#' @param plotParams Logical, whether to plot model parameters.
#' @param pltTitle Title for the plots.
#' @return A list with optimal lag summary and model outputs.
#' @export

analyzeLags <- function(f1="",f2="",dyad=c(),xname=f1,yname=f2, norm=F,window_size=60*5,window_step=window_size,start="", end="",func=fitDSModelForDyad,na.rm=T,simulate=F,dname=paste(xname,yname,sep="+"),measure="EDA", downsample=1, minLag=0,maxLag=5,noPlots=F,relativeToLag=-1,type=4,plotParams="raw",pltTitle=NULL) {
  lagList <- list()
n = 1;
  lags <- minLag:maxLag
  mdls <- list()
  for (lag in lags) {
    mdl <- analyzeDyad(f1,f2,dyad=dyad,xname=xname,yname=yname, norm=norm,window_size=window_size,window_step=window_step,start=start, end=end,func=func,na.rm=na.rm,simulate=simulate,dname=dname,lag=lag,noPlots = noPlots,plotParams = plotParams, pltTitle = pltTitle, measure = measure,downsample = downsample,type = type)
    mdlSummary <- mdl$summary
    mdlSummary$lag <- lag
    lagList[[n]] <- mdlSummary
    mdls[[n]] <- mdl
    n = n+1
  }
  
  output <- lagList[[1]]
  dxVals <- get.key(lagList,key = "dx.r.squared",as.list = T)
  dyVals <- get.key(lagList,key = "dy.r.squared",as.list = T)
  
  dxValsByLag <- unlist(do.call(cbind,dxVals))
  dyValsByLag <- unlist(do.call(cbind,dyVals))
  dxOptimalLag <- lags[apply(dxValsByLag,1,which.max)] #Compensate for lack of 0 indexing in R
  dyOptimalLag <- lags[apply(dyValsByLag,1,which.max)] #Compensate for lack of 0 indexing in R
  
  if(relativeToLag > -1){
    dxOptimal <- rowMaxs(dxValsByLag) - dxValsByLag[,relativeToLag+1]
    dyOptimal <- rowMaxs(dyValsByLag)- dyValsByLag[,relativeToLag+1]
  }
  else {
    dxOptimal <- rowMaxs(dxValsByLag)
    dyOptimal <- rowMaxs(dyValsByLag)
  }
  
  output$dx.r.squared <- dxOptimal
  output$dy.r.squared <- dyOptimal
  output$x.lag <- dxOptimalLag
  output$y.lag <- dyOptimalLag
  output$lag <- NULL
  
  print(plot.lagparams(output,xname=xname,yname=yname,title=paste(xname,"&",yname)))
  
  return(list(output=output,mdls=mdls))
  
}

#' Analyze dyad data for TVDSM analysis
#'
#' @description Processes two participant datasets or a dyad object, performs TVDSM analysis, generates plots, and returns model summaries.
#' @param f1 File path for participant 1.
#' @param f2 File path for participant 2.
#' @param dyad Optional dyad object.
#' @param xname Friendly name for participant 1.
#' @param yname Friendly name for participant 2.
#' @param norm Logical, whether to standardize data.
#' @param window_size Window size in seconds.
#' @param window_step Step size in seconds.
#' @param start Start time for analysis.
#' @param end End time for analysis.
#' @param func Function for model fitting (default fitDSModelForDyad).
#' @param na.rm Logical.
#' @param simulate Logical, whether to simulate a dyad.
#' @param dname Friendly name for the dyad.
#' @param lag Lag value for analysis.
#' @param noPlots Logical, if TRUE, no plots are produced.
#' @param plotParams Logical, whether to include model parameter plots.
#' @param pltTitle Title for the plots.
#' @param measure Measurement column name.
#' @param type Model type.
#' @param downsample Downsampling factor.
#' @param x.baseline Baseline for participant 1.
#' @param y.baseline Baseline for participant 2.
#' @param verbose Logical, whether to be verbose.
#' @param codes File path for condition codes.
#' @param useRealTime Logical, whether to interpret condition codes in real time.
#' @param incTSD Logical, whether to include time series descriptives.
#' @param ... Additional arguments.
#' @return A list containing TVDSM model outputs, summary, and plots.
#' @export

analyzeDyad <- function(f1="",f2="",dyad=c(),xname=f1,yname=f2, norm=F,window_size=60*5,window_step=window_size,start="", end="",func=fitDSModelForDyad,na.rm=T,simulate=F,dname=paste(xname,yname,sep="+"),lag=0,noPlots=F, plotParams="raw",pltTitle=paste(dname,"(","Lag","=",lag,")"), measure="EDA",type=4,downsample=1, x.baseline=NA, y.baseline=NA, verbose=F, codes="",useRealTime=F, incTSD=F, ...) {
  timeformat ="%Y-%m-%d %H:%M:%S"
  
  if(length(dyad) > 0){
    d <- dyad
  }
  else {
    p1 <- read.eda(f1)
    p2 <- read.eda(f2)
    
    if(simulate){
      d <- as.simulateddyad(p1,p2,norm=norm,cols = c(measure), verbose=verbose)
    }
    else {
      d <- as.dyad(p1,p2,norm=norm,cols = c(measure),verbose = verbose)
    }
    
    
  }
  
  if(is.na(x.baseline)) x.baseline <- mean(d[,2],na.rm=T)
  if(is.na(y.baseline)) y.baseline <- mean(d[,3],na.rm=T)
  fs <- getFS(d)
  
  FUN_TSD <- function(data){
    out <- func(data,lag=lag,type=type,downsample=downsample,x_mu=x.baseline,y_mu=y.baseline)
    out$xTSD <- tsDescriptives(data[,2],fs=fs)
    out$yTSD <- tsDescriptives(data[,3],fs=fs)
    return(out)
  }
  FUN <- function(data){
    return(func(data,lag=lag,type=type,downsample=downsample,x_mu=x.baseline,y_mu=y.baseline))
  }
  

  if(start != ""){
    start <- strptime(start,format=timeformat)
  }
  else {
    start <- min(d$Timestamp) 
  }
  if(end != ""){
    end <- strptime(end,format=timeformat)
  }
  else {
    end <- max(d$Timestamp)
  }
  d <- subset(d,((Timestamp >= start) & (Timestamp <= end)) )
  fs <- getFS(d)
  
  #Pre-downsample source data to match what is used in analysis
  ax <- decimate(d[,2], round(downsample*fs))
  ay <- decimate(d[,3],round(downsample*fs))
  t <- seq.int(1,to = length(d$Timestamp),by = round(downsample*fs))
  sourceData <- data.frame(Timestamp=t)
  sourceData[[names(d)[2]]] <- ax
  sourceData[[names(d)[3]]] <- ay
  
  
  if (window_size <= 0) {
    window_size <- length(d$Timestamp)/fs
    window_step <- window_size
  }
  if(incTSD){
    data <- o.window.list(d,window_size = window_size*fs, window_step=window_step*fs, FUN = FUN_TSD,na.rm = na.rm,verbose=verbose)
  }
  else {
    data <- o.window.list(d,window_size = window_size*fs, window_step=window_step*fs, FUN = FUN,na.rm = na.rm,verbose=verbose)
  }
  if(codes != ""){
    
    if (useRealTime) {
      ERCodes <- read.RTCodes(codes)
      validateConditionTimes(d,ERCodes) #Check that condition times don't exceed dyad start/end times and throw meaningful error messages
      
      start.ctime <- as.POSIXct(start,format=timeformat)
      end.ctime <- as.POSIXct(end,format=timeformat)
      if(verbose){print(paste("Dyad valid from ",start,"to",end,"with duration of",round(difftime(end.ctime,start.ctime,units="secs")),"seconds"))}
      ERCodes$Start <- difftime(as.POSIXct(ERCodes$Start.Time,format=timeformat),start.ctime, units = "secs")
      ERCodes$End <- difftime(as.POSIXct(ERCodes$End.Time,format=timeformat),start.ctime, units="secs")
    }
    
    else {
      ERCodes <- read.Codes(codes,...)
      ERCodes$Start  <- ERCodes$Start.Time
      ERCodes$End <- ERCodes$End.Time
    }
    l <- dim(ERCodes)[1]
    ERCodes[l+1,] <- c(NA,-1,-1,-1,-1)
    addNA(ERCodes$Condition)
    for(n in rev(seq_along(data))){
      if(is.na(data[[n]]) == F){
        dt <- data[[n]]$Start - start
        data[[n]]$Condition <- ERCodes$Condition[[l+1]]
        for (i in 1:l) {
          if( (dt >= ERCodes[i,"Start"]) & (dt < ERCodes[i,"End"])){
            data[[n]]$Condition <- ERCodes$Condition[[i]]
          }
        }
      }
    }
    
  }
  
  
  
  out <- list()
  
  if(noPlots == F){
    pltData<- plot.tvdsm(data,xname=xname,yname=yname,use.delta.rsquared = T,by.condition = (codes != ""),title = pltTitle,plotParams=plotParams,rawData=d, ...)
    print(pltData)
    out$plt <- pltData
    colnames(d) <- c("Timestamp",xname,yname)
    
    out$rawPlot <- plotDyad(d,ylabel = measure,title = pltTitle, ...)
    
  }
  n <- length(get.key(data,"x.r.squared"))
  mdls <- data
  mdlData <- data.frame(timestamp=as.POSIXct(get.key(mdls,"Timestamp"),origin = "1970-01-01"),
                        name=rep_len(dname,n),
                        Start=as.POSIXct(get.key(mdls,"Start"),origin = "1970-01-01"),
                        End=as.POSIXct(get.key(mdls,"End"),origin = "1970-01-01"),
                        dx.r.squared=get.key(mdls,"dx.r.squared"),
                        dy.r.squared=get.key(mdls,"dy.r.squared"),
                        x.selfreg=get.key(mdls,"b1"),
                        x.coreg=get.key(mdls,"b2"),
                        x.interaction=get.key(mdls,"b21"),
                        y.selfreg=get.key(mdls,"b4"),
                        y.coreg=get.key(mdls,"b5"),
                        y.interaction=get.key(mdls,"b45")
  )

  if(incTSD){
    df.xTSD <- ldply(mdls,.fun = function(d){return(as.data.frame(d$xTSD))})
    df.yTSD <- ldply(mdls,.fun = function(d){return(as.data.frame(d$yTSD))})
    colnames(df.xTSD) <- paste0("x.",colnames(df.xTSD))
    colnames(df.yTSD) <- paste0("y.",colnames(df.yTSD))
    mdlData <- as.data.frame(cbind(mdlData,df.xTSD,df.yTSD))
  }
  out$mdls <- mdls
  out$summary <- mdlData
  out$sourceData <- sourceData
  return(out)
}

#' Describe dyad interdependence
#'
#' @description Computes descriptive statistics for dyad interdependence metrics based on TVDSM outputs.
#' @param mdls A list of TVDSM model outputs.
#' @param xname Friendly name for participant 1.
#' @param yname Friendly name for participant 2.
#' @return A summary of interdependence measures.
#' @export

describeDyad <- function(mdls,xname="Participant 1", yname="Participant 2") {
  dx <- get.key(mdls, "dx.r.squared")
  dy <- get.key(mdls, "dy.r.squared")
  dyad <- rowMeans(cbind(dx,dy))
  
  group <- c(rep_len(xname,length(dx)),rep_len(yname,length(dy)),rep_len("Dyad",length(dy)))
  interdependence <- c(dx,dy,dyad)
  data <- data.frame(Participant=group,I3=interdependence)
  
  return(o.describeBy(interdependence, group,tex = F,digits = 3))
}

#' Save interpersonal analysis data to a CSV file
#'
#' @description Exports the summary data (and optionally condition-based data) from a dyad analysis to a CSV file.
#' @param dyadData A list containing TVDSM dyad analysis output.
#' @param outputFilename Optional filename for the output CSV.
#' @param I3.only Logical, if TRUE, only saves the I3 metric.
#' @param aname Friendly name for participant 1.
#' @param bname Friendly name for participant 2.
#' @param by.condition Logical, if TRUE, includes condition information.
#' @return The data frame that was saved.
#' @export

saveInterpersonalData <- function(dyadData, outputFilename=NULL, I3.only=F,aname=NULL,bname=NULL,by.condition=F){
  timeformat ="%Y-%m-%dT%H:%M:%S%z"
  d <- dyadData$mdls
  tz <- attr(d[[1]]$Timestamp,"tz")
  if( (is.null(d[[1]]$Condition) == F) & by.condition ){
    if(I3.only){
      data <- data.frame(Timestamp=strftime(as.POSIXlt(as.numeric(get.key(d, "Timestamp" )),origin = "1970-01-01",tz = tz), format=timeformat),Condition=get.key(d,"Condition"), a.r.squared=get.key(d,"dx.r.squared"),b.r.squared=get.key(d,"dy.r.squared"))
      
    }
    else {
      data <- data.frame(Timestamp=strftime(as.POSIXlt(as.numeric(get.key(d, "Timestamp" )),origin = "1970-01-01",tz = tz), format=timeformat),
                         Condition=get.key(d,"Condition"))
      
      data[paste0(aname,".I3")]=get.key(d,"dx.r.squared")
      data[paste0(bname,".I3")]=get.key(d,"dy.r.squared")
      data[paste0(aname,".SelfReg")]=get.key(d,"b1")
      data[paste0(aname,".CoReg")]=get.key(d,"b2")
      data[paste0(bname,".SelfReg")]=get.key(d,"b4")
      data[paste0(bname,".CoReg")]=get.key(d,"b5")
    }
  }
  else {
    if(I3.only){
      a.name <- paste0("P-",aname,"_I3")
      b.name <- paste0("P-",bname,"_I3")
      data <- data.frame(Timestamp=strftime(as.POSIXlt(as.numeric(get.key(d, "Timestamp" )),origin = "1970-01-01",tz = tz), format=timeformat))
      data[a.name] <- get.key(d,"dx.r.squared")
      data[b.name] <- get.key(d,"dy.r.squared")
    }
    else {
      data <- data.frame(Timestamp=strftime(as.POSIXlt(as.numeric(get.key(d, "Timestamp" )),origin = "1970-01-01",tz = tz), format=timeformat))
      data[paste0(aname,".I3")]=get.key(d,"dx.r.squared")
      data[paste0(bname,".I3")]=get.key(d,"dy.r.squared")
      data[paste0(aname,".SelfReg")]=get.key(d,"b1")
      data[paste0(aname,".CoReg")]=get.key(d,"b2")
      data[paste0(bname,".SelfReg")]=get.key(d,"b4")
      data[paste0(bname,".CoReg")]=get.key(d,"b5")
    }
    
  }
  
  if(is.null(outputFilename)){
    outputFilename <- paste(as.character(dyadData$summary$name[[1]]),".csv",sep="")
  }
  abnames <- strsplit(as.character(dyadData$summary$name[[1]]),split = "+")
  #data$a.name <- abnames[1]
  #data$b.name <- abnames[2]

  
  write.csv(data,file = outputFilename,row.names=FALSE,quote = F)
  return(data)
}

#' Save data for graphing from TVDSM analysis
#'
#' @description Exports the R-squared summary metrics of the TVDSM analysis to CSV files.
#' @param mdls A list containing TVDSM model outputs and summary data.
#' @param xfname Filename for saving participant X's R-squared data.
#' @param yfname Filename for saving participant Y's R-squared data.
#' @return None.
#' @export

saveDataForGraphing <- function(mdls,xfname="Dyad_X.csv",yfname="Dyad_X.csv"){
  timeformat ="%Y-%m-%d %H:%M:%S"
  data <- mdls$summary
  data$Timestamp <- strftime(as.POSIXlt(as.numeric(data$timestamp),origin = "1970-01-01"), format=timeformat)
  write.csv(data.frame(Timestamp=data$Timestamp,R2=data$dx.r.squared),file = xfname)
  write.csv(data.frame(Timestamp=data$Timestamp, R2=data$dy.r.squared),file = yfname)
}

#' Compute dynamical correlation using bootstrap methods
#'
#' @description Calculates the dynamical correlation between participant data using a bootstrap approach.
#' @param a_files A vector of file paths for participant A.
#' @param b_files A vector of file paths for participant B.
#' @param read_func Function to read data files (default read.eda).
#' @param cols A vector of measurement columns (default "EDA").
#' @param sr Sampling rate (default 32).
#' @param randPairs Logical, if TRUE, uses random pairing.
#' @return Bootstrapped dynamical correlation test results.
#' @export

dynamicalCorrelation <- function(a_files, b_files,read_func=read.eda,cols=c("EDA"),sr=32,randPairs=F){
  n <- length(a_files)
  dyads <- list()
  nrows <- 0
  
  for (i in 1:n) {
    a.raw <- read_func(a_files[[i]])
    b.raw <- read_func(b_files[[i]])
    d <- as.dyad(a.raw,b.raw,cols=cols)
    dyads[[i]] <- d
    if(nrows < nrow(d)){nrows = nrow(d)}
  }
  a_mat <- matrix(data=NA, nrow = nrows, ncol = n)
  b_mat <- matrix(data=NA, nrow = nrows, ncol = n)
  for(j in 1:n){
    d <-dyads[[j]] 
    a_mat[1:nrow(d),j] <-d[1:nrow(d),2]
    b_mat[1:nrow(d),j] <-d[1:nrow(d),3]
  }
  t <-seq(from=0,to=nrows/sr,length.out=nrows)
  
  #results.test <- ind_DC(a_mat,b_mat,t,na=T)
  results.boot <- boot_test_DC(a_mat,b_mat,t,ms = T,randPairs = randPairs)
  
  return(results.boot)
}

#' Generate an I3 matrix for interdependence analysis
#'
#' @description Constructs a 3D matrix of I3 values across participants and time windows and optionally saves it as JSON.
#' @param data A nested list of dyad analysis results.
#' @param outputFile Optional filename to save the JSON output.
#' @return A list containing the legend, I3 matrix, and timestamps.
#' @export

generateI3Mat <- function(data, outputFile=NULL){
  timestamps = data[[1]][[1]]$summary$timestamp
  nWindows <- length(data[[1]][[1]]$summary$dx.r.squared)
  ids <- unique(union(names(data),unlist(lapply(data, names))))
  
  nParticipants <- length(ids)
  
  mat <- array(data = 0,dim = c(nParticipants,nParticipants,nWindows),dimnames = list(P1=ids,P2=ids,Time=timestamps))
  tsd <- list()
  for (P1 in names(data)) {
    P1.data <- data[[P1]]
    for (P2 in names(P1.data)){
      P2.data <- P1.data[[P2]]
      print(paste0("mat[",P1,",",P2,",]"))
      #Swap dx and dy so that the I3 values for a given row represent effect of P_row on P_col
      mat[P1,P2,] <- P2.data$summary$dy.r.squared
      mat[P2,P1,] <- P2.data$summary$dx.r.squared
      if(dim(P2.data$summary)[2] > 8){
        tsd[[P1]] <- P2.data$summary[,13:23]
        tsd[[P2]] <- P2.data$summary[,24:34]
        colnames(tsd[[P1]]) <- str_sub(colnames(tsd[[P1]]),start=3)
        colnames(tsd[[P2]]) <- str_sub(colnames(tsd[[P2]]),start=3)
      }
    }
  }
  
  mat[,,][mat < 0] <- 0
  output = list(legend=ids,I3Mat=mat,timestamps=timestamps)
  if(length(tsd) > 0){
    output$tsd <- tsd
  }

  if(is.null(outputFile) == F){
    write(toJSON(output,matrix = "columnmajor",pretty = T,na = NULL,factor = "string"),file = outputFile)
  }
  
  return(output)
}

#' Analyze a group of participant data files
#'
#' @description Performs pairwise TVDSM analysis on a group of participant data files.
#' @param files A vector of file paths or a directory path containing CSV files.
#' @param window_size Window size in seconds for TVDSM analysis.
#' @param window_step Window step size in seconds.
#' @param readFunction Function to read each data file (default read.eda).
#' @param scanDirs Logical, if TRUE, treats files as directory and scans for CSV files.
#' @param ... Additional arguments passed to the analysis functions.
#' @return A list containing analysis results for each participant pair.
#' @export

analyzeGroup <- function(files,window_size=60,window_step=10,readFunction=read.eda,scanDirs=T,...){
  if(dir.exists(files) & scanDirs){
    files <-list.files(path = files,pattern = "\\.csv$",full.names = T)
  }
  files <- sort(files,decreasing = T)
  filenames <- basename(files)
  nTotal = 0
  nFiles = length(files) 
  output <- list()
  ids <- filenames
  physioData <- list()
  for (f in files) {
    physioData[[f]] <- readFunction(f)
  }
  startTime <- as.POSIXct(max(unlist(lapply(physioData, function(d){return(min(d$Timestamp))}))),origin="1970-01-01")
  endTime <- as.POSIXct(min(unlist(lapply(physioData, function(d){return(max(d$Timestamp))}))),origin="1970-01-01")
  print(paste0("Overlap Start: ",strftime(startTime,format = "%Y-%m-%d %H:%M:%S")," | Overlap End: ",strftime(endTime,format = "%Y-%m-%d %H:%M:%S")))
  for(i in 1:(nFiles-1)) {
    personA = files[i]
    aName = ids[i]
    idx2 = i+1
    output[[aName]] <- list()
    for(j in idx2:nFiles){
      personB = files[j]
      bName = ids[j]
      nTotal = nTotal+1
      d <- as.dyad(physioData[[personA]],physioData[[personB]])
      d <- subset(d,((Timestamp >= startTime) & (Timestamp <= endTime)) )
      if(j > i){
        try(output[[aName]][[bName]]<- analyzeDyad(dyad = d,xname=aName,yname=bName,window_size = window_size,window_step = window_step, type = 4,lag = 0,plotParams = "raw",measure="EDA",...))
      }
    }
  }
  
  return(output)
}

#' Combine participant data into a group dyad
#'
#' @description Merges overlapping time series data from multiple participants into a single data frame.
#' @param pn A list of participant data frames.
#' @param cols A vector of column names to include (default "EDA").
#' @param norm Logical, whether to normalize data.
#' @param na.interpolate Logical, if TRUE interpolates missing values.
#' @param interpolation.method Function used for interpolation (default na.spline).
#' @return A combined data frame of group time series data.
#' @export

as.group <- function(pn,cols=c("EDA"),norm=F,na.interpolate=T,interpolation.method=na.spline) {
  pn.name <- names(pn)
  pn.start <- max(unlist(lapply(pn, FUN=function(p){return(min(p$Timestamp))})))
  pn.end <- min(unlist(lapply(pn,FUN=function(p){return(max(p$Timestamp))})))
  
  if(pn.end <= pn.start){
    print("Error: no overlap found between timeseries")
    return()
  }
  print(paste("Dyad Overlap: ",pn.start,"to",pn.end))
  pn.cropped <- lapply(pn, FUN=function(p){return( subset(p,((Timestamp >= pn.start) & (Timestamp < pn.end)) ) )})
  dyad <- data.frame(Timestamp=pn.cropped[[1]]$Timestamp)
  
  for(i in 1:length(pn)){
    p <- pn.cropped[[i]]
    
    for(c in cols) {
      cname <- paste(pn.name[[i]],c,sep = ".")
      
      if(norm){
        dyad[[cname]] <- o.scale(p[[c]])
      }
      else {
        dyad[[cname]] <- p[[c]]
      }
      
      if(na.interpolate){
        dyad[[cname]] <- interpolation.method(dyad[[cname]])
        
      }
      
    }
  }  
  return(dyad)
  
}

o.Lag <- function(x, lag=0, robust=T) {
  if (!is.numeric(x)) {
    warning("Input vector is not numeric. Converting to numeric.")
    x <- as.numeric(x)
  }
  
  if(robust){
    dx <- (Lag(x,(-1*lag)-1) - Lag(x,(-1*lag)+1))/2
  }
  else {
    dx <- Lag(x,(-1*lag)-1) -  Lag(x,(-1*lag))
  }
  
  return(dx)
}


