library(lattice)
library(lme4)
library(psychometric)
library(reshape)
library(e1071)
library(lmerTest)
library(cowplot)
library(lmtest)

library(systemfit)

library(progress)


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

#' Perform TVDSM analysis on a dyad by condition
#'
#' @param f1 path to csv file with data for participant 1
#' @param f2 path to csv file with data for participant 2
#' @param codes path to csv file with condition timing data
#' @param useRealTime treat condition times listed in condition timing sheet as calendar times, e.g. "1/2/16 13:23:45"
#' @param type the type of TVDSM analysis to use
#' @param cols a vector indicating which columns to select from each file when creating a dyad
#' @param dname a string indicating a friendly name for the dyad (used in plots and summary data frame)
#' @param p1.name a string indicating a friendly name for the participant whose data was loaded at path indicated by \code{f1}
#' @param p2.name a string indicating a friendly name for the participant whose data was loaded at path indicated by \code{f2}
#' @param plotParams should model coefficients be ploted
#' @param lag the lag to be used in the TVDSM analysis (in seconds)
#' @param downsample the desired sample period in seconds after downsampling (e.g. 8Hz = 0.125, 2Hz = 0.5, 1hz= 1.0)
#' @param lag the lag to be used in the TVDSM analysis (in seconds)
#' @return a list containing models (mdls), plot (plt), and a summary data frame (summary)
#' @export
#' @examples
#' analyzeByCondition("PersonA.csv","PersonB.csv",codes="ConditionTimes.csv")


analyzeByCondition <- function(f1="",f2="", dyad=c(),codes="",useRealTime=F, type=4,cols=c("EDA"),dname="Dyad",p1.name="Participant 1",p2.name="Participant 2",lag=0,plotParams=T,downsample=1,func=computeStateSpace, verbose=F,start="", end="",  timeformat ="%Y-%m-%d %H:%M:%S"){
  
  if(length(dyad) > 0){
    d <- dyad
  }
  else {
    p1 <- read.eda(f1)
    p2 <- read.eda(f2)
    d <- as.dyad(p1,p2,norm = T,cols=cols)
    
  }

  

  if(start != ""){
    start <- strptime(start,format=timeformat)
  }
  else {
    start <- d$Timestamp[[1]]
  }
  if(end != ""){
    end <- strptime(end,format=timeformat)
  }
  else {
    end <- d$Timestamp[[length(d$Timestamp)]]
  }
  d <- subset(d,((Timestamp >= start) & (Timestamp < end)) )
 
  
  if (useRealTime) {
    ERCodes <- read.RTCodes(codes)
    validateConditionTimes(d,ERCodes) #Check that condition times don't exceed dyad start/end times and throw meaningful error messages
    
    start.ctime <- as.POSIXct(start,format=timeformat)
    end.ctime <- as.POSIXct(end,format=timeformat)
    if(verbose){print(paste("Dyad valid from ",start,"to",end,"with duration of",round(difftime(end.ctime,start.ctime,units="secs")),"seconds"))}
    ERCodes$Start.Time <- difftime(as.POSIXct(ERCodes$Start.Time,format=timeformat),start.ctime, units = "secs")
    ERCodes$End.Time <- difftime(as.POSIXct(ERCodes$End.Time,format=timeformat),start.ctime, units="secs")
  }
  else {
    ERCodes <- read.Codes(codes)
  }
  
  if(verbose){View(ERCodes)}

  FS <- getFS(d)
  mdls <- list()
  
  
  x.baseline <- mean(d[,2],na.rm=T)
  y.baseline <- mean(d[,3],na.rm=T)
  
  FUN <- function(data){
    return(func(data,lag=lag,type=type,downsample=downsample,x_mu=x.baseline,y_mu=y.baseline,verbose=verbose));
  }
  if(verbose){
    pb <- progress_bar$new(total = dim(ERCodes)[[1]])
  }
  for(n in 1:dim(ERCodes)[[1]]){
    if(!is.na(ERCodes$Condition[[n]])){
      start = ERCodes$Start.Time[[n]]*FS
      end = ERCodes$End.Time[[n]]*FS
      if(end >= dim(d)[[1]]){
        end = dim(d)[[1]]-1
      }
      mdls[[n]] <- FUN(d[start:end,])
      mdls[[n]]$Condition <- ERCodes$Condition[[n]]
      mdls[[n]]$Duration <- (ERCodes$End.Time[[n]] - ERCodes$Start.Time[[n]])
      mdls[[n]]$Start <- ERCodes$Start.Time[[n]] 
      mdls[[n]]$End <- ERCodes$End.Time[[n]] 
      mdls[[n]]$Timestamp <- ERCodes$Start.Time[[n]] + (ERCodes$End.Time[[n]] - ERCodes$Start.Time[[n]])/2
      mdls[[n]]$dyad.name <- dname 
      
    }
    if(verbose){
      pb$tick()
    }
    
  }
  startD = ERCodes$Start.Time[[1]]*FS
  endD = ERCodes$End.Time[[dim(ERCodes)[1]]]*FS
  
  rawD <- d[start:end,]
  
  plt <- plot.ssparams(mdls,xname = p1.name,yname=p2.name,use.delta.rsquared = T, title=dname,plotParams = plotParams,rawData = rawD)
  print(plt)
  mdlData <- data.frame(timestamp=get.key(mdls,"Timestamp"),name=rep_len(dname,n),Condition=get.key(mdls,"Condition"), Start=get.key(mdls,"Start"),End=get.key(mdls,"End"), x.r.squared=get.key(mdls,"dx.r.squared"),y.r.squared=get.key(mdls,"dy.r.squared"),x.selfreg=get.key(mdls,"b1"),x.coreg=get.key(mdls,"b2"),x.interaction=get.key(mdls,"b21"),y.selfreg=get.key(mdls,"b4"),y.coreg=get.key(mdls,"b5"),y.interaction=get.key(mdls,"b45"))
  
  output <- list(mdls=mdls,summary=mdlData,plt=plt)
  
  return(output)
  
}




get.key <- function(lst, key,as.list=F) {
  if(as.list == T){
    return(as.list(lapply(lst,FUN = function(l){return(try(l[[key]]))}),recursive = F ))
  }
  else{
    return(unlist( lapply(lst,FUN = function(l){return(try(l[[key]]))}),recursive = F ))
  }
}

plot.lagparams <- function(data,xname="x",yname="y",title=paste(xname,"&",yname),save=F,f="plot.pdf",use.delta.rsquared=T,by.condition=T, autoscale=T) {
  data <- na.omit(data)
  r2.labels <- c(bquote(Delta~R[.(xname)]^2),bquote(Delta~R[.(yname)]^2))
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


plot.ssparams <- function(data,xname="x",yname="y",title=paste(xname,"&",yname),save=F,f="plot.pdf",use.delta.rsquared=T,by.condition=T, autoscale=T,plotParams=T,rawData=NULL) {
  data <- na.omit(data)
  r2.labels <- c(bquote(R[.(xname)]^2),bquote(R[.(yname)]^2))
  if(use.delta.rsquared){
    x.r2 <- get.key(data, key="dx.r.squared")
    y.r2 <- get.key(data, key="dy.r.squared")
    r2.labels <- c(bquote(Delta~R[.(xname)]^2),bquote(Delta~R[.(yname)]^2))
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
  
  if(is.na(params.df$b21) & is.na(params.df$b45)){
    param.labels <- c(bquote(.(xname)[S-R]), bquote(.(xname)[Co-R]),bquote(.(yname)[S-R]),bquote(.(yname)[Co-R]))
    glegend1 <- scale_colour_discrete(name  = "Beta Coefficents",
                                      breaks=c("b1","b2", "b4","b5"),
                                      labels=param.labels)
    
  }
  else {
    param.labels <- c(bquote(.(xname)[S-R]), bquote(.(xname)[Co-R]), bquote(.(xname)[I]), bquote(.(yname)[S-R]),bquote(.(yname)[Co-R]),bquote(.(yname)[I]))
    
  glegend1 <- scale_colour_discrete(name  = "Coefficients",
                                      breaks=c("b1","b2","b21", "b4","b5","b45"),
                                      labels=param.labels)
  }

  glegend2 <- scale_colour_discrete(name  = "Effect Size",
                                    breaks=c("x.r2","y.r2"),
                                    labels=r2.labels)
  gstyle <- guides(colour = guide_legend(
    label.theme = element_text(size=15,angle=0),label.hjust=0,label.vjust=0),title.theme=element_text(size=15,angle=0))
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
  plt1 <- ggplot(data=params1, aes(x=Timestamps,y=value,colour=variable)) + geom_line(size=1) + xlab(x_label) + ylab(NULL) + ggtitle(title.1) +glegend1 +gstyle 
  
  plt2 <- ggplot(data=params2, aes(x=Timestamps,y=value,colour=variable)) + geom_line(size=1) + xlab(x_label) + ylab(NULL) + ggtitle(title.2) + glegend2 + gstyle
  
  if(point_size > 1.25) {
    plt1 <- plt1 + geom_point(size=point_size) 
    plt2 <- plt2 + geom_point(size=point_size) 
  }
  
  if(by.condition){
    plt2 <- plt2 + geom_rect(aes(fill=Condition,ymin=-0.05,ymax=0.0,xmin=start,xmax=end),linetype=0,alpha=0.5) 
  }
  plt1 <- plt1 + scale_x_datetime(limits=c(x_min, x_max))
  plt2 <- plt2 + scale_x_datetime(limits=c(x_min, x_max)) 
  if(plotParams == T){
    combinedPlt <- plot_grid(plt1,plt2,align = "hv",ncol = 1,nrow = 2,rel_heights = c(2,1))
  }
  else if(plotParams == "raw"){
    plt.raw <- plotDyad(rawData,title = title,ylabel = "")
    combinedPlt <- plot_grid(plt.raw,plt2,align = "hv",ncol = 1,nrow = 2,rel_heights = c(2,1))
  }
  else {
    combinedPlt <- plt2
  }

  if(save){
    ggsave(f)
  }
  return(combinedPlt)
  
}



plot.density <- function(data,group="Dyad", xname='Particpant 1 Beta', yname="Participant 2 Beta") {
  x.beta <- get.key(data, key="b2")
  y.beta <- get.key(data, key="b4")

  kp1 <- kde2d(x.beta,y.beta,n=500)
  contour(kp1, col='blue',)
  title(paste("Kernel Density Plots -", group), ylab=yname, xlab=xname)
  
  gg <- ggplot(data=kp1,mapping = aes(x,y,z)) + geom_contour()
  print(gg)
}


statespace.fiml <- function(x,y,x_mu=mean(x,na.rm=T),y_mu=mean(y,na.rm=T),type=2, step=0.25,p.value=0.01,verbose=F,lag=0){
  if(is.null(lag)){
    x_prime <- Lag(x,-1*lag)
    y_prime <- Lag(y,-1*lag)
    
  }
  else {
    x_prime <-o.Lag(x,lag=lag)
    y_prime <- o.Lag(y,lag=lag)
    
  }

  s1 <- ""
  s2 <- ""
  
  data <- data.frame(x_prime=x_prime,y_prime=y_prime, x=x,y=y)
  
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
    
    x_model <- x_prime ~ I(x_mu-x) + y_prime
    y_model <- y_prime ~ I(y_mu-y) + x_prime
    
    
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


computeStateSpace <- function(dyad,type=4,downsample=1,lag=0,x_mu=NULL,y_mu=NULL,verbose=F) {
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
  
  tryCatch(mdl <- statespace.fiml(ax,ay,p.value = 0.05,type=type,lag=lag*newFS,x_mu = x_mu,y_mu=y_mu,verbose=verbose),
           error=function(e){
             print(e)
           })
  mdl$Timestamp <- min(dyad$Timestamp)
  mdl$Duration <- (max(dyad$Timestamp) - min(dyad$Timestamp))
  mdl$Start <- min(dyad$Timestamp)
  mdl$End <- max(dyad$Timestamp)

  return(mdl)
}

urlForModel <- function(m) {
  a1 <- m$b1
  a2 <- m$b2
  b1 <- m$b4
  b2 <- m$b5
  minV <- min(m$mat[,1:2],na.rm = T)
  maxV <- max(m$mat[,1:2], na.rm = T)
  return(paste("file:///Users/OliverWS/git/DCS%20Visualization/index.html?a1=", a1,"&a2=",a2,"&b1=",b1,"&b2=",b2,"&min=",minV,"&max=",maxV,sep=""))
}


analyzeLags <- function(f1="",f2="",dyad=c(),xname=f1,yname=f2, norm=F,window_size=60*5,window_step=window_size,start="", end="",func=computeStateSpace,na.rm=T,simulate=F,dname=paste(xname,yname,sep="+"),measure="EDA", downsample=1, minLag=0,maxLag=5,noPlots=F,relativeToLag=-1,type=4,plotParams=T,pltTitle=NULL) {
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



#' Perform TVDSM analysis on a dyad by condition
#'
#' @param f1 path to csv file with data for participant 1
#' @param f2 path to csv file with data for participant 2
#' @param dyad if a valid dyad object is provided (produced via \code{\link{as.dyad}}), analysis will be run on the dyad directly, and \code{f1} and \code{f2} will be ignored
#' @param type the type of TVDSM analysis to use
#' @param dname a string indicating a friendly name for the dyad (used in plots and summary data frame)
#' @param xname a string indicating a friendly name for the participant whose data was loaded at path indicated by \code{f1}
#' @param yname a string indicating a friendly name for the participant whose data was loaded at path indicated by \code{f2}
#' @param norm should participant data be standardized prior to TVDSM analysis
#' @param window_size the size of the sliding window to use for TVDSM analysis, in seconds
#' @param window_step the increment to be used when sliding \code{window_size}, in seconds
#' @param measure a string indicating which column to select from each file when creating a dyad
#' @param pltTitle pass a string to override plot title
#' @param plotParams should model coefficients be ploted
#' @param noPlots if true, plots will not be generated
#' @param lag the lag to be used in the TVDSM analysis (in seconds)
#' @param downsample the desired sample period in seconds after downsampling (e.g. 8Hz = 0.125, 2Hz = 0.5, 1hz= 1.0)
#' @param x.baseline provide a value to be treated as baseline for participant 1 in TVDSM analysis. Defaults to mean of participant 1's data.
#' @param y.baseline provide a value to be treated as baseline for participant 2 in TVDSM analysis. Defaults to mean of participant 2's data.
#' @param downsample the desired sample period in seconds after downsampling (e.g. 8Hz = 0.125, 2Hz = 0.5, 1hz= 1.0)
#' @return a list containing models (mdls), plot (plt), and a summary data frame (summary)
#' @export
#' @examples
#' #Perform TVDSM analysis using a 5 minute sliding window incremented every 10 second, with a lag of 5 seconds
#' analyzeDyad("PersonA.csv","PersonB.csv",window_size=5*60, window_step=10,lag=5)



analyzeDyad <- function(f1="",f2="",dyad=c(), xname=f1,yname=f2, norm=F,window_size=60*5,window_step=window_size,start="", end="",func=computeStateSpace,na.rm=T,simulate=F,dname=paste(xname,yname,sep="+"),lag=0,noPlots=F, plotParams=T,pltTitle=paste(dname,"(","Lag","=",lag,")"), measure="EDA",type=4,downsample=1, x.baseline=NA, y.baseline=NA, verbose=F, codes="",useRealTime=F, incTSD=F) {
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
  d <- subset(d,((Timestamp >= start) & (Timestamp < end)) )
  
  fs <- getFS(d)
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
      ERCodes$Start.Time <- difftime(as.POSIXct(ERCodes$Start.Time,format=timeformat),start.ctime, units = "secs")
      ERCodes$End.Time <- difftime(as.POSIXct(ERCodes$End.Time,format=timeformat),start.ctime, units="secs")
    }
    
    else {
      ERCodes <- read.Codes(codes)
    }
    l <- dim(ERCodes)[1]
    ERCodes[l+1,] <- c(NA,-1,-1)
    addNA(ERCodes$Condition)
    for(n in 1:length(data)){
      if(is.na(data[[n]]) == F){
        dt <- data[[n]]$Start - start
        data[[n]]$Condition <- ERCodes$Condition[[l+1]]
        for (i in 1:l) {
          if( (dt >= ERCodes[i,"Start.Time"]) & (dt < ERCodes[i,"End.Time"])){
            data[[n]]$Condition <- ERCodes$Condition[[i]]
          }
        }
      }
    }
    
  }
  
  
  
  out <- list()
  
  if(noPlots == F){
    pltData<- plot.ssparams(data,xname=xname,yname=yname,use.delta.rsquared = T,by.condition = (codes != ""),title = pltTitle,plotParams=plotParams,rawData=d)
    print(pltData)
    out$plt <- pltData
    out$rawPlot <- plotDyad(d,ylabel = measure,title = pltTitle)
    
  }
  n <- length(get.key(data,"x.r.squared"))
  mdls <- data
  mdlData <- data.frame(timestamp=get.key(mdls,"Timestamp"),start=get.key(mdls,"Start"),end=get.key(mdls,"End"),name=rep_len(dname,n), x.r.squared=get.key(mdls,"x.r.squared"),y.r.squared=get.key(mdls,"y.r.squared"),dx.r.squared=get.key(mdls,"dx.r.squared"),dy.r.squared=get.key(mdls,"dy.r.squared"))
  out$mdls <- mdls
  out$summary <- mdlData
  return(out)
}

describeDyad <- function(mdls,xname="Participant 1", yname="Participant 2") {
  dx <- get.key(mdls, "dx.r.squared")
  dy <- get.key(mdls, "dy.r.squared")
  dyad <- rowMeans(cbind(dx,dy))
  
  group <- c(rep_len(xname,length(dx)),rep_len(yname,length(dy)),rep_len("Dyad",length(dy)))
  interdependence <- c(dx,dy,dyad)
  data <- data.frame(Participant=group,R2=interdependence)
  
  return(o.describeBy(interdependence, group,tex = F,digits = 3))
}

saveInterpersonalData <- function(dyadData, outputFilename=NULL, I3.only=F,aname="a",bname="b",by.condition=F){
  timeformat ="%Y-%m-%dT%H:%M:%S%z"
  d <- dyadData$mdls
  tz <- attr(d[[1]]$Timestamp,"tz")
  if( (is.null(d[[1]]$Condition) == F) & by.condition ){
    if(I3.only){
      data <- data.frame(Timestamp=strftime(as.POSIXlt(as.numeric(get.key(d, "Timestamp" )),origin = "1970-01-01",tz = tz), format=timeformat),Condition=get.key(d,"Condition"), a.r.squared=get.key(d,"dx.r.squared"),b.r.squared=get.key(d,"dy.r.squared"))
      
    }
    else {
      data <- data.frame(Timestamp=strftime(as.POSIXlt(as.numeric(get.key(d, "Timestamp" )),origin = "1970-01-01",tz = tz), format=timeformat),Condition=get.key(d,"Condition"), a.r.squared=get.key(d,"dx.r.squared"),b.r.squared=get.key(d,"dy.r.squared"),a.selfreg=get.key(d,"b1"),a.coreg=get.key(d,"b2"),a.interaction=get.key(d,"b21"),b.selfreg=get.key(d,"b4"),b.coreg=get.key(d,"b5"),b.interaction=get.key(d,"b45"))
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
      data <- data.frame(Timestamp=strftime(as.POSIXlt(as.numeric(get.key(d, "Timestamp" )),origin = "1970-01-01",tz = tz), format=timeformat), a_I3=get.key(d,"dx.r.squared"),b_I3=get.key(d,"dy.r.squared"),a_selfreg=get.key(d,"b1"),a_coreg=get.key(d,"b2"),a_interaction=get.key(d,"b21"),b_selfreg=get.key(d,"b4"),b_coreg=get.key(d,"b5"),b_interaction=get.key(d,"b45"))
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

saveDataForGraphing <- function(mdls,xfname="Dyad_X.csv",yfname="Dyad_X.csv"){
  timeformat ="%Y-%m-%d %H:%M:%S"
  data <- mdls$summary
  data$Timestamp <- strftime(as.POSIXlt(as.numeric(data$timestamp),origin = "1970-01-01"), format=timeformat)
  write.csv(data.frame(Timestamp=data$Timestamp,R2=data$dx.r.squared),file = xfname)
  write.csv(data.frame(Timestamp=data$Timestamp, R2=data$dy.r.squared),file = yfname)
}



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