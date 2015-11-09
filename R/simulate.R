library(Hmisc)
library(readxl)
library(lattice)
library(lme4)
library(psychometric)
library(reshape)
library(coefplot)
library(e1071)
library(lmerTest)
CBSL.update <- function(){
  library(devtools)
  install_github("OliverWS/CBSL.R",auth_token = "2f83b9ce082240ab6a07075eaee6ff32e1c954f8")
}
simulateDyad <- function(duration=600,fs=32,selfReg.coef=0.5,coReg.coef=1.0,interaction.coef=0,lag=0,mu=1,sd=2,sr.ratio=0.5) {
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
  timestamps <- seq(from = start, by=dt , length.out = duration*fs)
  
  x.t <- seq_len(length(timestamps))
  x.noise <- o.scale(rnorm(n=length(x.t), mean=mu, sd=0.1))
  x.signal <-o.scale(o.smooth(o.scale(rnorm(n=length(x.t), mean=mu, sd=1)),window = fs*30))
  x <- x.signal*sr.ratio #+ x.noise*(1.0-sr.ratio)

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
  lines(y.t,y,col="green")
  outputData <- data.frame(Timestamp = timestamps)
  outputData[,2] <- x
  outputData[,3] <- y
  colnames(outputData) <- c("Timestamp","Person A","Person B")
  outputData <- na.omit(outputData)
  data <- melt(as.data.frame(outputData),id.vars = c("Timestamp"))
  g1<- ggplot(data = data,mapping = aes(x=Timestamp,y=value,col=variable)) + geom_line(size=1) + xlab("Time") + ylab("EDA") + ggtitle(eq)
  g1 <- g1 + scale_y_continuous()
  g1 <- g1 + scale_x_datetime()
  print(g1)
  
  return(outputData)
}

simulatePartner <- function(x.signal,fs=32,selfReg.coef=0.5,coReg.coef=1.0,interaction.coef=0,lag=0,mu=1,sd=2,sr.ratio=0.5) {
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
  lines(y.t,y,col="green")
  outputData <- data.frame(Timestamp = timestamps)
  outputData[,2] <- x
  outputData[,3] <- y
  colnames(outputData) <- c("Timestamp","Person A","Person B")
  outputData <- na.omit(outputData)
  data <- melt(as.data.frame(outputData),id.vars = c("Timestamp"))
  g1<- ggplot(data = data,mapping = aes(x=Timestamp,y=value,col=variable)) + geom_line(size=1) + xlab("Time") + ylab("EDA") + ggtitle(eq)
  g1 <- g1 + scale_y_continuous()
  g1 <- g1 + scale_x_datetime()
  print(g1)
  
  return(outputData)
}



bootstrapSimulation <- function(n=100) {
  r<-replicate(n,computeStateSpace(simulateDyad(duration=600,fs = 1,sr.ratio =0.95,lag=0,mu = 0.5,sr = 0,cr =0.25,i = 0),downsample = 1,x_mu = 0.5,y_mu=0.5,lag = 0,type = 2),simplify = "array")
  return(unlist(r[11,]))
}

examples <- function(){
  setwd("~/Dropbox/Dynamical Systems Modeling Paper/")
  s <- simulateDyad(duration = 300,fs = 1.0,mu = 0.5,sr.ratio = 0.95)
  x <- s[,2]
  sr <- 0
  i <- 0
  n = 1
  for(cr in seq(0,0.25,0.05)){
    sp <- simulatePartner(x.signal = x,fs = 1.0,mu = 0.5,sr.ratio = 1.0,cr = cr,sr=sr,i = 0)
    ggsave(paste("Frame_",n,".png",sep=""),width = 10,height = 8)
    n = n + 1
    #ggsave(paste("Simulation_SR=",sr,"_CR=",cr,"_I=",i,".png",sep=""))
  }
  for(sr in seq(0,0.25,0.05)){
    sp <- simulatePartner(x.signal = x,fs = 1.0,mu = 0.5,sr.ratio = 1.0,cr = cr,sr=sr,i = 0)
    ggsave(paste("Frame_",n,".png",sep=""),width = 10,height = 8)
    n = n + 1
    #ggsave(paste("Simulation_SR=",sr,"_CR=",cr,"_I=",i,".png",sep=""))
  }
}




