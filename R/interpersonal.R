source("~/git/OWS.R/OWS.R")
source("~/git/CBSL.R/R/utils.R")
source("~/git/CBSL.R/R/read.R")
source("~/git/CBSL.R/R/plot.R")
source("~/git/CBSL.R/R/physio.R")
library(lattice)
library(lme4)
library(psychometric)
library(reshape)
library(e1071)
library(lmerTest)
library(cowplot)


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



analyzeByCondition <- function(f1="",f2="",codes="",type=2,cols=c("EDA"),dname="Dyad.1",p1.name="Infant",p2.name="Parent",lag=0){
  
  ERCodes <- read.Codes(codes)
  View(ERCodes)
  
  
  p1 <- read.eda(f1)
  p2 <- read.eda(f2)
  
  d <- as.dyad(p1,p2,norm = T,cols=cols)
  plot(d$Timestamp,d[,2],type="l",col="blue") 
  lines(d$Timestamp,d[,3],col="red")
  FS <- getFS(d)
  mdls <- list()
  for(n in 1:dim(ERCodes)[[1]]){
    if(!is.na(ERCodes$Condition[[n]])){
      start = ERCodes$Start.Time[[n]]*FS
      end = ERCodes$End.Time[[n]]*FS
      if(end >= dim(d)[[1]]){
        end = dim(d)[[1]]-1
      }
      print(paste("Start",start,"End",end))
      mdls[[n]] <- computeStateSpace(d[start:end,],type = type,downsample =1,lag=lag )
      mdls[[n]]$Condition <- ERCodes$Condition[[n]]
      mdls[[n]]$Duration <- (ERCodes$End.Time[[n]] - ERCodes$Start.Time[[n]])
      mdls[[n]]$Start <- ERCodes$Start.Time[[n]] 
      mdls[[n]]$End <- ERCodes$End.Time[[n]] 
      mdls[[n]]$Timestamp <- ERCodes$Start.Time[[n]] + (ERCodes$End.Time[[n]] - ERCodes$Start.Time[[n]])/2
      mdls[[n]]$dyad.name <- dname 
      
    }
    
    
  }
  

  
  plot.ssparams(mdls,xname = p1.name,yname=p2.name,use.delta.rsquared = T, title=dname)
  
  
  
  
  return(mdls)
  
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


plot.ssparams <- function(data,xname="x",yname="y",title=paste(xname,"&",yname),save=F,f="plot.pdf",use.delta.rsquared=T,by.condition=T, autoscale=T) {
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
  
  params1 <- melt(subset(params.df,select=c("Timestamps","b1","b2", "b4","b5","b21","b45")),id.vars = c("Timestamps"))
  params2 <- melt(subset(params.df,select=c("Timestamps","x.r2","y.r2","Condition","start","end")),id.vars = c("Timestamps","start","end","Condition"))
  
  glegend1 <- scale_colour_discrete(name  ="Variable",
                                    breaks=c("b1","b2","b21", "b4","b5","b45"),
                                    labels=c(expression(alpha[S-R]), expression(alpha[Co-R]),expression(alpha[I]),expression(beta[S-R]),expression(beta[Co-R]),expression(beta[I])))
  
  glegend2 <- scale_colour_discrete(name  ="Variable",
                                    breaks=c("x.r2","y.r2"),
                                    labels=r2.labels)
  gstyle <- guides(colour = guide_legend(
    label.theme = element_text(size=15,angle=0),label.hjust=0,label.vjust=0),title.theme=element_text(size=15,angle=0))
  x_min = min(start)
  x_max = max(end)
  plt1 <- ggplot(data=params1, aes(x=Timestamps,y=value,colour=variable)) + geom_line(size=1) + geom_point(size=3) + xlab("Time") + ylab("Coefficient") + ggtitle(title) +glegend1 +gstyle 
  plt2 <- ggplot(data=params2, aes(x=Timestamps,y=value,colour=variable)) + geom_line(size=1) + geom_point(size=3) + xlab("Time") + ylab(bquote(R^2)) + ggtitle(title) + glegend2 + gstyle
  if(by.condition){
    plt2 <- plt2 + geom_rect(aes(fill=Condition,ymin=-0.05,ymax=0.0,xmin=start,xmax=end),linetype=0,alpha=0.5)
  }
  plt1 <- plt1 + scale_x_datetime(limits=c(x_min, x_max))
  plt2 <- plt2 + scale_x_datetime(limits=c(x_min, x_max)) 
  combinedPlt <- plot_grid(plt1,plt2,align = "hv",ncol = 1,nrow = 2,rel_heights = c(2,1))
  
  
  
  if(save){
    ggsave(f)
  }
  return(combinedPlt)
  
}

computeStateSpace <- function(dyad,type=2,downsample=1,lag=0) {
  FS <- getFS(dyad)
  ax <- decimate(dyad[,2], downsample*FS)
  ay <- decimate(dyad[,3],downsample*FS)
  mdl <- statespace.fiml(ax,ay,p.value = 0.05,type=type,lag=lag)
  mdl$Timestamp <- dyad[1,"Timestamp"]
  mdl$Duration <- (max(mdl$Timestamp) - min(mdl$Timestamp))
  mdl$Start <- mdl$Timestamp[[1]] 
  mdl$End <- max(mdl$Timestamp)
  
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


analyzeLags <- function(f1,f2,xname=f1,yname=f2, norm=F,window_size=60*5,window_overlap=0,start="", end="",func=computeStateSpace,na.rm=T,simulate=F,dname=paste(xname,yname,sep="+"),minLag=0,maxLag=5,noPlots=F,relativeToLag=-1) {
  lagList <- list()
  n = 1;
  lags <- minLag:maxLag
  for (lag in lags) {
    mdl <- analyzeDyad(f1,f2,xname=xname,yname=yname, norm=norm,window_size=window_size,window_overlap=window_overlap,start=start, end=end,func=func,na.rm=na.rm,simulate=simulate,dname=dname,lag=lag,noPlots = noPlots)
    mdlSummary <- mdl$summary
    mdlSummary$lag <- lag
    lagList[[lag+1]] <- mdlSummary
    n = n+1
  }
  
  output <- lagList[[minLag+1]]
  dxVals <- get.key(lagList,key = "dx.r.squared",as.list = T)
  dyVals <- get.key(lagList,key = "dy.r.squared",as.list = T)
  
  dxValsByLag <- unlist(do.call(cbind,dxVals))
  dyValsByLag <- unlist(do.call(cbind,dyVals))
  dxOptimalLag <- apply(dxValsByLag,1,which.max) -1 #Compensate for lack of 0 indexing in R
  dyOptimalLag <- apply(dyValsByLag,1,which.max) -1 #Compensate for lack of 0 indexing in R
  
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
  
  return(output)
  
}


analyzeDyad <- function(f1,f2,xname=f1,yname=f2, norm=F,window_size=60*5,window_overlap=0,start="", end="",func=computeStateSpace,na.rm=T,simulate=F,dname=paste(xname,yname,sep="+"),lag=0,noPlots=F) {
  if(simulate){
    p1 <- f1
    p2 <- f2
  }
  else {
    p1 <- read.eda(f1)
    p2 <- read.eda(f2)
  }
  timeformat ="%Y-%m-%d %H:%M:%S"
  FUN <- function(data){
    return(func(data,lag=lag));
  }
  
  
#  if(simulate){
#    d <- as.simulateddyad(p1,p2,norm=norm)
#  }
#  else {
    d <- as.dyad(p1,p2,norm=norm)
#  }
  
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
  
  data <- o.window.list(d,window_size = window_size*fs, window_overlap=window_overlap*fs, FUN = FUN,na.rm = na.rm)
  out <- list()
  
  if(noPlots == F){
    pltData<- plot.ssparams(data,xname=xname,yname=yname,use.delta.rsquared = T,by.condition = F,title = paste(xname,"+",yname,"(","Lag","=",lag,")"))
    pltData
    out$plt <- pltData
  }
  
  mdls <- data
  n <- length(get.key(mdls,"x.r.squared"))
  mdlData <- data.frame(timestamp=get.key(mdls,"Timestamp"),name=rep_len(dname,n), x.r.squared=get.key(mdls,"x.r.squared"),y.r.squared=get.key(mdls,"y.r.squared"),dx.r.squared=get.key(mdls,"dx.r.squared"),dy.r.squared=get.key(mdls,"dy.r.squared"))
  out$mdls <- data
  out$summary <- mdlData
  return(out)
}


saveData <- function(mdls,xfname="Dyad_X.csv",yfname="Dyad_X.csv"){
  timeformat ="%Y-%m-%d %H:%M:%S"
  data <- mdls$summary
  data$Timestamp <- strftime(as.POSIXlt(as.numeric(data$timestamp),origin = "1970-01-01"), format=timeformat)
  write.csv(data.frame(Timestamp=data$Timestamp,R2=data$dx.r.squared),file = xfname)
  write.csv(data.frame(Timestamp=data$Timestamp, R2=data$dy.r.squared),file = yfname)
}

