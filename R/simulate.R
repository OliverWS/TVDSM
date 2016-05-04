library(Hmisc)
library(matrixStats)
library(signal)
library(readxl)
library(lattice)
library(lme4)
library(psychometric)
library(reshape)
library(coefplot)
library(e1071)
library(lmerTest)
library(progress)

simulateDyad <- function(duration=600,fs=32,selfReg.coef=0.5,coReg.coef=1.0,interaction.coef=0,lag=0,mu=1,sd=2,trend=0,sr.ratio=0.95) {
  sr = selfReg.coef
  cr = coReg.coef
  i = interaction.coef
  
  ds.predict <- function(x,y,type=2) {
    if(type == 2){
      dx = 0.005 + sr*(mu-x) + cr*(y-x) + i*(mu-x)*(y-x)
    }
    return(dx)
  }
  
  start=strptime("2015-01-01 0:0:0","%Y-%m-%d %H:%M:%S")
  
  dt <- as.difftime((1.0/fs),units = "secs")
  timestamps <- seq(from = start, by=dt , length.out = duration*fs)
  
  x.t <- seq_len(length(timestamps))
  x.noise <- o.scale(rnorm(n=length(x.t), mean=mu, sd=0.1))
  x.signal <-o.scale(o.smooth(o.scale(rnorm(n=length(x.t), mean=mu, sd=1)),window = fs*30))
  x.signal <- x.signal + seq(0.0,trend,length.out=length(x.signal))
  x <- x.signal*sr.ratio #+ x.noise*(1.0-sr.ratio)

  y.t <- x.t
  y.signal <- rep_len(x.signal[1],length.out = length(x.t))
  
  for(n in 1:(length(x.t)-1)){
    y.signal[n+1] = y.signal[n]+ds.predict(y.signal[n],x.signal[n])
  }
  eq <- bquote(Model: frac(Delta~x,Delta~t) == .(sr)%.%(bar(x) - x[t]) + .(cr)%.%(y[t]-x[t]) + .(i)%.%(bar(x) - x[t]) %*% (y[t]-x) )
  y.signal <- Lag(y.signal,shift = -lag)

  y.noise <-  o.scale(rnorm(n=length(x.t), mean=mu, sd=0.1))
  y <- y.signal*sr.ratio + y.noise*(1.0-sr.ratio) + runif(n = 1)
  print(length(timestamps))
  print(length(x))
  print(length(y))
  outputData <- data.frame(Timestamp = timestamps)
  outputData[,2] <- x
  outputData[,3] <- y
  colnames(outputData) <- c("Timestamp","Person A","Person B")
  outputData <- na.omit(outputData)
  data <- melt(as.data.frame(outputData),id.vars = c("Timestamp"))
  g1<- ggplot(data = data,mapping = aes(x=Timestamp,y=value,col=variable)) + geom_line(size=1) + xlab("Time") + ylab("EDA") # + ggtitle(eq)
  g1 <- g1 + scale_y_continuous()
  g1 <- g1 + scale_x_datetime()

  return(list(plt=g1,data=outputData,eq=eq))
}

simulatePartner <- function(x.signal,fs=32,selfReg.coef=0.5,coReg.coef=1.0,interaction.coef=0,lag=0,mu=1,sd=2,sr.ratio=0.95) {
  sr = selfReg.coef
  cr = coReg.coef
  i = interaction.coef
  
  
  ds.predict <- function(x,y,type=2) {
    if(type == 2){
      dx = sr*(mu-x) + cr*(y-x) + i*(mu-x)*(y-x)
    }
    return(dx)
  }
  
  start=strptime("2015-01-01 0:0:0","%Y-%m-%d %H:%M:%S")
  
  dt <- as.difftime((1.0/fs),units = "secs")
  timestamps <- seq(from = start, by=dt , length.out = length(x.signal))
  
  x.t <- seq_len(length(timestamps))
  x <- x.signal
  y.t <- x.t
  y.signal <- rep_len(x.signal[1],length.out = length(x.t))
  
  for(n in 1:(length(x.t)-1)){
    y.signal[n+1] = y.signal[n]+ds.predict(y.signal[n],x.signal[n])
  }
  eq <- bquote(Model: frac(Delta~x,Delta~t) == .(sr)%.%(bar(x) - x[t]) + .(cr)%.%(y[t]-x[t]) + .(i)%.%(bar(x) - x[t]) %*% (y[t]-x) )
  y.signal <- Lag(y.signal,shift = -lag)
  
  y.noise <-  o.scale(rnorm(n=length(x.t), mean=mu, sd=0.1))
  y <- y.signal*sr.ratio + y.noise*(1.0-sr.ratio)
  print(length(timestamps))
  print(length(x))
  print(length(y))
  outputData <- data.frame(Timestamp = timestamps)
  outputData[,2] <- x
  outputData[,3] <- y
  colnames(outputData) <- c("Timestamp","Person A","Person B")
  outputData <- na.omit(outputData)
  data <- melt(as.data.frame(outputData),id.vars = c("Timestamp"))
  g1<- ggplot(data = data,mapping = aes(x=Timestamp,y=value,col=variable)) + geom_line(size=1) + xlab("Time") + ylab("EDA") #+ ggtitle(eq)
  g1 <- g1 + scale_y_continuous()
  g1 <- g1 + scale_x_datetime()
  print(g1)

  return(list(data=outputData,plt=g1,eq=eq))
}
o.scale <- function(x,use.sd=F){
  if(use.sd){
    sx <- (x - min(x, na.rm=T))/sd(x,na.rm = T)
    
  }
  else {
    sx <- (x - min(x, na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))    
  }
  return(sx)
  
}

o.smooth <- function(x,window=10){
  mat <- matrix(nrow=length(x),ncol = window)
  for(n in 1:window){
    l <- n - window/2
    mat[,n] <- Lag(x,l)
  }
  return(rowWeightedMeans(mat,w=gausswin(window),na.rm = T))
  
}


bootstrapSimulation <- function(n=100) {
  r<-replicate(n,computeStateSpace(simulateDyad(duration=600,fs = 1,sr.ratio =0.95,lag=0,mu = 0.5,sr = 0,cr =0.25,i = 0),downsample = 1,x_mu = 0.5,y_mu=0.5,lag = 0,type = 2),simplify = "array")
  return(unlist(r[11,]))
}

sampleSignal <- function(){
  data <- read_excel(paste("~/Dropbox/CouplesData/All Raw Data/","002",".xlsx",sep=""))
  
  d.sc <- makeDyad(cbind(data$`SC-A_`,data$`SC-C_`),sample_rate = 0.1)
  
  return(d.sc[,2])
}

makeDyad <- function(edaData, sample_rate=0.1, start=strptime("2015-01-01 0:0:0","%Y-%m-%d %H:%M:%S"),norm=T) {
  dt <- as.difftime(sample_rate,units = "secs")
  timestamps <- seq(from = start, by=dt , length.out = dim(edaData)[1])
  outputData <- data.frame(Timestamp = timestamps)
  if(norm) {
    outputData[,2] <- o.scale(edaData[,1])
    outputData[,3] <- o.scale(edaData[,2])
    
  }
  else {
    outputData[,2] <- edaData[,1]
    outputData[,3] <- edaData[,2]
  }
  return(na.omit(outputData))
}

compareSim <- function(x.signal = sampleSignal(), fs=0.1,coReg.coef=0.1,selfReg.coef=0.05){
  sp <- simulatePartner(x.signal = x.signal,fs =0.1,mu = 0.5,sr.ratio = 1.0,coReg.coef = coReg.coef,selfReg.coef=selfReg.coef,i = 0)
  a <- analyzeDyad(dyad=sp$data,xname = "Person A",yname="Person B",window_size = 600)
  plot_grid(sp$plt,a$plt,ncol=1,nrow=2)
}

examples <- function(){
  setwd("~/Dropbox/Dynamical Systems Modeling Paper/")
  data <- read_excel(paste("~/Dropbox/CouplesData/All Raw Data/","002",".xlsx",sep=""))
  
  d.sc <- makeDyad(cbind(data$`SC-A_`,data$`SC-C_`),sample_rate = 0.1)
  
  x <- d.sc[,2]
  sr <- 0
  i <- 0
  n = 1
  for(cr in seq(0,0.25,0.05)){
    sp <- simulatePartner(x.signal = x,fs = 1.0,mu = 0.5,sr.ratio = 1.0,coReg.coef = cr,selfReg.coef=sr,i = 0)
    #ggsave(paste("Frame_",n,".png",sep=""),width = 10,height = 8)
    n = n + 1
    #ggsave(paste("Simulation_SR=",sr,"_CR=",cr,"_I=",i,".png",sep=""))
  }
  for(sr in seq(0,0.25,0.05)){
    sp <- simulatePartner(x.signal = x,fs = 1.0,mu = 0.5,sr.ratio = 1.0,cr = cr,sr=sr,i = 0)
    #ggsave(paste("Frame_",n,".png",sep=""),width = 10,height = 8)
    n = n + 1
    #ggsave(paste("Simulation_SR=",sr,"_CR=",cr,"_I=",i,".png",sep=""))
  }
}


runFalsePairings <- function(participantAList, participantBList, window_size=60, window_step=60, downsample=1, lag=0) {
  n <- length(participantAList)
  nFalse = 1
  truePairs <- list()
  falsePairs <- list()
  nPossiblePairs = n*n*1.0
  nTotal = 0
  
  pb <- progress_bar$new(total = nPossiblePairs)
  
  
  for(i in 1:n) {
    personA = participantAList[[i]]
    aName = personA
    personB = participantBList[[i]]
    bName = personB
    #print(paste("True:",aName,"<=>",bName))
    
    truePairs[[i]] <- analyzeDyad(f1 = personA, f2=personB,xname = aName, yname = bName, window_size = window_size, window_step = window_step, lag = lag,noPlots = T)
    for(j in 1:n){
      if(j != i){
        personB = participantBList[[j]]
        bName = personB
        nTotal = nTotal+1
        if(j > i){
          #print(paste("False:",aName,"<=>",bName))
          falsePairs[[nTotal]] <- analyzeDyad(f1 = personA, f2=personB,xname = aName, yname = bName, simulate = T, window_size = window_size, window_step = window_step, lag = lag)
        }
      }
      pb$tick()
    }
    
  }
  
  trueDataX <- lapply(get.key(get.key(truePairs, "summary",as.list = T), "dx.r.squared", as.list = T), mean)
  trueDataY <- lapply(get.key(get.key(truePairs, "summary",as.list = T), "dy.r.squared",as.list = T), mean)
  falseDataX <- lapply(get.key(get.key(falsePairs, "summary",as.list = T), "dx.r.squared", as.list = T), mean)
  falseDataY <- lapply(get.key(get.key(falsePairs, "summary",as.list = T), "dy.r.squared",as.list = T), mean)
  trueData <- c(unlist(trueDataX), unlist(trueDataY))
  falseData <- c(unlist(falseDataX), unlist(falseDataY))
  print(t.test(trueData,falseData,))
  close(pb)
  return(list(truePairs=truePairs, falsePairs=falsePairs))
  
}
  
  