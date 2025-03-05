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
library(stringr)
library(compute.es)

#' Simulate a Dyadic Time Series
#' 
#' This function simulates a dyadic time series for two interacting persons using self- and co-regulation parameters.
#' 
#' @param duration Duration of simulation in seconds (default: 600).
#' @param fs Sampling frequency (default: 32).
#' @param selfReg.coef Self-regulation coefficient (default: 0.5).
#' @param coReg.coef Co-regulation coefficient (default: 1.0).
#' @param interaction.coef Interaction coefficient (default: 0).
#' @param lag Lag to apply (default: 0).
#' @param mu Mean value for baseline signal (default: 1).
#' @param sd Standard deviation for noise (default: 2).
#' @param trend Linear trend added to the signal (default: 0).
#' @param sr.ratio Ratio of signal to noise (default: 0.95).
#' @return A list with elements: plt (ggplot object), data (simulated data frame), and eq (model equation).
#' @export
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

#' Simulate a Partner Time Series
#'
#' This function simulates a partner's time series based on an input signal and regulation parameters.
#'
#' @param x.signal Numeric vector representing the input signal.
#' @param fs Sampling frequency (default: 32).
#' @param selfReg.coef Self-regulation coefficient (default: 0.5).
#' @param coReg.coef Co-regulation coefficient (default: 1.0).
#' @param interaction.coef Interaction coefficient (default: 0).
#' @param lag Lag to apply (default: 0).
#' @param mu Mean value (default: 1).
#' @param sd Standard deviation (default: 2).
#' @param trend Linear trend to add (default: 0).
#' @param sr.ratio Ratio of signal to noise (default: 0.95).
#' @return A list with elements: data (simulated partner data frame), plt (ggplot object), and eq (model equation).
#' @export
simulatePartner <- function(x.signal,fs=32,selfReg.coef=0.5,coReg.coef=1.0,interaction.coef=0,lag=0,mu=1,sd=2,trend=0,sr.ratio=0.95) {
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

#' Bootstrap Dyadic Simulation
#'
#' Runs a bootstrapping simulation using the dyadic simulation and model estimation.
#'
#' @param n Number of bootstrap replications (default: 100).
#' @return A vector of bootstrapped simulation statistics.
#' @export
bootstrapSimulation <- function(n=100) {
  r<-replicate(n,fitDSModelForDyad(simulateDyad(duration=600,fs = 1,sr.ratio =0.95,lag=0,mu = 0.5,sr = 0,cr =0.25,i = 0),downsample = 1,x_mu = 0.5,y_mu=0.5,lag = 0,type = 2),simplify = "array")
  return(unlist(r[11,]))
}

#' Sample Signal
#'
#' Reads a sample signal from an Excel file for simulation purposes.
#'
#' @return A numeric vector representing the sample signal.
#' @export
sampleSignal <- function(){
  data <- read_excel(paste("~/Dropbox/CouplesData/All Raw Data/","002",".xlsx",sep=""))
  
  d.sc <- makeDyad(cbind(data$`SC-A_`,data$`SC-C_`),sample_rate = 0.1)
  
  return(d.sc[,2])
}

#' Create Dyad Data Frame
#'
#' Constructs a dyadic data frame from EDA data with generated timestamps.
#'
#' @param edaData Matrix or data frame with two columns representing signals.
#' @param sample_rate Sampling rate (default: 0.1).
#' @param start Start time as POSIXlt (default: 2015-01-01 0:0:0).
#' @param norm Logical indicating whether to normalize data (default: TRUE).
#' @return A data frame with columns: Timestamp, Person A, and Person B.
#' @export
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

#' Compare Simulated Signals
#'
#' Compares simulated signals using partner simulation and dyadic analysis, then plots the results.
#'
#' @param x.signal Signal for person A (default: sampleSignal()).
#' @param fs Sampling frequency (default: 0.1).
#' @param coReg.coef Co-regulation coefficient (default: 0.1).
#' @param selfReg.coef Self-regulation coefficient (default: 0.05).
#' @return A plot grid comparing simulation outputs.
#' @export
compareSim <- function(x.signal = sampleSignal(), fs=0.1,coReg.coef=0.1,selfReg.coef=0.05){
  sp <- simulatePartner(x.signal = x.signal,fs =0.1,mu = 0.5,sr.ratio = 1.0,coReg.coef = coReg.coef,selfReg.coef=selfReg.coef,i = 0)
  a <- analyzeDyad(dyad=sp$data,xname = "Person A",yname="Person B",window_size = 600)
  plot_grid(sp$plt,a$plt,ncol=1,nrow=2)
}

#' Run Example Simulations
#'
#' Executes example simulation scenarios over a range of regulation coefficients.
#' @export
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

#' Run False Pairings Analysis
#'
#' Runs TVDSM analysis on true and false pairings from lists of participant signals.
#'
#' @param participantAList List of signals for person A.
#' @param participantBList List of signals for person B.
#' @param window_size Window size for analysis (default: 60).
#' @param window_step Step size between windows (default: window_size).
#' @param downsample Downsample factor (default: 1).
#' @param lag Lag to be applied (default: 0).
#' @param type Analysis type specifier (default: 3).
#' @return A list containing results for true pairs and false pairs.
#' @export
runFalsePairings <- function(participantAList, participantBList, window_size=60, window_step=window_size, downsample=1, lag=0,type=3) {
  n <- length(participantAList)
  nFalse = 1
  truePairs <- list()
  falsePairs <- list()
  nPossiblePairs = n*n
  nTotal = 0
  print(paste("Running TVDSM on",n,"true pairs and",nPossiblePairs,"false pairs",sep = " "))
  pb <- progress_bar$new(
    format = "Processing :what :spin [:bar] :percent eta: :eta",
    clear = F, total = nPossiblePairs)
  msg <- ""
  pb$tick(0)
  for(i in 1:n) {
    personA = participantAList[[i]]
    aName = personA
    for(j in 1:n){
      personB = participantBList[[j]]
      bName = personB
      msg <- paste(aName,"<=>",bName)
      pb$tick(tokens = list(what=msg))
      
      if(j != i){
        nTotal = nTotal+1
        if(j > i){
          falsePairs[[nTotal]] <- analyzeDyad(f1 = personA, f2=personB,xname = aName, yname = bName, simulate = T, window_size = window_size, window_step = window_step, lag = lag,type = type)
        }
      }
      else {
        truePairs[[i]] <- analyzeDyad(f1 = personA, f2=personB,xname = aName, yname = bName, window_size = window_size, window_step = window_step, lag = lag,noPlots = F,type = type)
      }
    }
    
  }
  out <- list(truePairs=truePairs, falsePairs=falsePairs)
  compareSimulatedDyads(out)
  return(out)
  
}

#' Compare Simulated Dyads
#'
#' Compares true dyad pairings with simulated (false) dyad pairings using statistical tests and visualization.
#'
#' @param pairs A list with elements 'truePairs' and 'falsePairs'.
#' @param title Title for the comparison plot (default: "True vs. Simulated Dyads").
#' @return A list containing the compiled data frame, plot object, and t-test results.
#' @export
compareSimulatedDyads <- function(pairs,title="True vs. Simulated Dyads"){
  truePairs <- pairs$truePairs
  falsePairs <- pairs$falsePairs
  customMean <- function(x){
    return(mean(x,na.rm=T))
  }
  xname <- function(a){
    try(return(str_split(a$name[[1]],pattern = "(\\+)")[[1]][1]))
    return(NA)
  }
  yname <- function(a){
    try(return(str_split(a$name[[1]],pattern = "(\\+)")[[1]][2]))
    return(NA)
  }
  
  trueDataX <- lapply(get.key(get.key(truePairs, "summary",as.list = T), "dx.r.squared", as.list = T), customMean)
  trueDataXName <- lapply(get.key(truePairs, "summary",as.list = T), xname)
  trueDataY <- lapply(get.key(get.key(truePairs, "summary",as.list = T), "dy.r.squared",as.list = T), customMean)
  trueDataYName <- lapply(get.key(truePairs, "summary",as.list = T), yname)
  
  falseDataX <- lapply(get.key(get.key(falsePairs, "summary",as.list = T), "dx.r.squared", as.list = T), customMean)
  falseDataXName <- lapply(get.key(falsePairs, "summary",as.list = T), xname)
  
  falseDataY <- lapply(get.key(get.key(falsePairs, "summary",as.list = T), "dy.r.squared",as.list = T), customMean)
  falseDataYName <- lapply(get.key(falsePairs, "summary",as.list = T), yname)
  
  
  #trueData <- na.omit(rowMeans(cbind(unlist(trueDataX), unlist(trueDataY)),na.rm = T))
  #falseData <- na.omit(rowMeans(cbind(unlist(falseDataX), unlist(falseDataY)),na.rm = T))

  trueData <- c(unlist(trueDataX),unlist(trueDataY))
  falseData <- c(unlist(falseDataX),unlist(falseDataY)) 
  
  
  
  
  df <- rbind.data.frame(data.frame(value=trueData, type="True Pair", name=c(unlist(trueDataXName),unlist(trueDataYName))),data.frame(value=falseData,type="Simulated Pair", name=c(unlist(falseDataXName),unlist(falseDataYName))))
  df <- na.omit(df)
  
  print(sd(trueData,na.rm = T))
  print(sd(falseData,na.rm = T))
  t <- t.test(trueData,falseData,var.equal = T)
  print(t)
  es = compute.es::tes(t$statistic,length(trueData),length(falseData))
  t$es <- es
  print(paste0("Cohen's d = ",es$d))
  print(paste("Total pairs:",nrow(df)))
  
  o.hist(df$value,group = df$type)
  
  
  p <- ggplot(df, aes(type, value,col=type)) + geom_point()
  p +   stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", width = 0.3) + ggtitle(title)
  return(list(df=df,plt=p,ttest=t))
}
  