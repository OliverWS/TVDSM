require(signal)
require(matrixStats)
require(Hmisc)
require(lavaan)

#' Compute Root Mean Square of Successive Differences (RMSSD)
#'
#' @param hr.period A numeric vector of heart rate period values.
#' @return The RMSSD value.
#' @export
RMSSD <- function(hr.period) {
  hr.period.delta <- diff(hr.period)
  return(sqrt(mean(hr.period.delta*hr.period.delta)))
}

#' Compute Linear Trend Slope
#'
#' @param data A numeric vector.
#' @return The slope of the linear trend.
#' @export
linearTrend <- function(data) {
  t <- seq_along(data)
  mdl <- lm(data ~ t)
  trend <- coef(mdl)[[2]]
  return(trend)
}

#' Compute Standard Error of a Vector
#'
#' @param x A numeric vector.
#' @return The standard error.
#' @export
se <- function(x) sqrt(var(x)/length(x))

#' Compute Time Series Descriptive Statistics
#'
#' @param data Numeric vector of data.
#' @param fs Sample rate (default is 32).
#' @return A list of descriptive statistics including mean, se, sd, var, trend, min, max, slope, slopeRatio, range and autocorrelation.
#' @export
tsDescriptives <- function(data,fs=32) {
  m <- mean(data, na.rm=T)
  st.err <- se(data)
  st.dev <- sd(data,na.rm=T)
  #Compute linear trend
  trend<-NA
  try((trend <- linearTrend(data)*fs*60),silent = T)
  max.data <- max(data)
  min.data <- min(data)
  range.data <- range(data)
  slope <- sum(diff(data))
  slopeRatio <- sum(diff(data) > 0)/length(data)
  try(data.ds <- decimate(data,q = fs),silent=T)
  ac <- NA
  try(ac <- acf(data.ds,lag.max=1,type="correlation",plot = F)$acf[2], silent=T)
  
  return(list(mean=m,se=st.err, sd=st.dev, var=var(data,na.rm=T), trend=trend,min=min.data,max=max.data,slope=slope,slopeRatio=slopeRatio,range=(max.data-min.data),autocorrelation=ac))
}

#' Dummy Code By Condition
#'
#' @param filename Filename for the EDA data.
#' @param codes Either a path to a file with codes or a dataframe with codes.
#' @param outputFile Output file name (default is filename appended with "_DummyCoded.csv").
#' @param saveFile Logical indicating whether to save the dummy coded file.
#' @param codeFunc Function to read codes, default is read.RTCodes.
#' @param startCol Column name for start time, default "Start.Time".
#' @param endCol Column name for end time, default "End.Time".
#' @param conditionCol Column name for condition, default "Condition".
#' @param timeformat Time format for conversion.
#' @return Data frame with dummy coded data.
#' @export
dummyCodeByCondition <- function(filename,codes,outputFile=paste(filename,"_DummyCoded.csv",sep=""),saveFile=F,codeFunc=read.RTCodes,startCol="Start.Time",endCol="End.Time",conditionCol="Condition",timeformat="%Y-%m-%d %H:%M:%S") {
  
  if(typeof(codes) == "character"){
    conditions <- codeFunc(codes)
  }
  else {
    conditions <- codes
  }
  conditions[[startCol]] <-  strptime(conditions[[startCol]], format=timeformat)
  conditions[[endCol]] <-  strptime(conditions[[endCol]], format=timeformat)
  eda <- read.eda(file = filename)
  fs <- getFS(eda)
  nConditions <- dim(conditions)[1]
  tz <- attr(eda$Timestamp[1],"tz")
  eda$Condition <- NA
  for (i in 1:nConditions) {
    start <- conditions[[startCol]][i]
    end <- conditions[[endCol]][i]
    condition <- as.character(conditions[[conditionCol]][i])
    validTimes <- ((eda$Timestamp >= start) & (eda$Timestamp < end))
    eda$Condition[which(validTimes)] <- condition
  }
  if(saveFile){
    # for (i in 2:nConditions) {
    #   condition.2 <- as.character(conditions[["Condition"]][i])
    #   condition.1 <- as.character(conditions[["Condition"]][i-1])
    #   sub <- subset(eda, ((eda$Condition == condition.2) | (eda$Condition == condition.1)))
    #   write.csv(sub,file = paste(filename,"_DummyCoded_",i-1,"&",i,".csv"),quote = T)
    # 
    # }
    eda$Timestamp <- strftime(eda$Timestamp, format="%Y-%m-%dT%H:%M:%OS3")
    write.csv(as.data.frame(eda), file = outputFile,col.names = T,row.names = F,)

  }
  return(eda)
}

#' Descriptives By Condition
#'
#' @param filename Path to the EDA data file.
#' @param codes Codes data (either file path or dataframe).
#' @param col Column for analysis (default "EDA").
#' @param title Title for the plot.
#' @param codeFunc Function to read codes, default is read.RTCodes.
#' @param startCol Start time column name, default "Start.Time".
#' @param endCol End time column name, default "End.Time".
#' @param conditionCol Condition column name, default "Condition".
#' @param timeformat Time format string.
#' @return Data frame of descriptive statistics.
#' @export
descriptivesByCondition <- function(filename, codes,col="EDA",title=filename,codeFunc=read.RTCodes,startCol="Start.Time",endCol="End.Time",conditionCol="Condition",timeformat="%Y-%m-%d %H:%M:%S") {
  
  if(typeof(codes) == "character"){
    conditions <- codeFunc(codes)
  }
  else {
    conditions <- codes
  }
  eda <- read.eda(file = filename)
  conditions[[startCol]] <-  strptime(conditions[[startCol]], format=timeformat)
  conditions[[endCol]] <-  strptime(conditions[[endCol]], format=timeformat)
  
  fs <- getFS(eda)
  nConditions <- dim(conditions)[1]
  results <- list()
  tz <- attr(eda$Timestamp[1],"tz")
  for (i in 1:nConditions) {
    start <- conditions[[startCol]][i]
    end <- conditions[[endCol]][i]
    condition <- conditions[[conditionCol]][i]
    validTimes <- ((eda$Timestamp >= start) & (eda$Timestamp < end))
    data <- subset(eda,validTimes)
    desc <- tsDescriptives(data[[col]])
    
    results[[i]] <- desc
  }

  results <- as.data.frame(t(as.data.frame(sapply(results, unlist))))
  rownames(results) <- NULL
  #Fix weird issue where values get turned into text
  results <- lapply(results,as.double)
  results$Condition <- conditions[[conditionCol]]
  results$Start <- conditions[[startCol]]
  results$End <- conditions[[endCol]]
  
  #g2 <-ggplot(data = subset(eda, ((eda$Timestamp >= min(results$Start)) & (eda$Timestamp < max(results$End)))) ,mapping = aes(x=Timestamp,y=EDA)) + geom_line(size=1,col="#1FBFC4") + xlab("Time") + ylab("EDA")

  try(g1 <- plotTSDescriptives(results,title=title))
  
  #plt <- plot_grid(g1,g2,ncol=1,nrow=2,align = "h")
  try(print(g1))
  return(results)
  
  
}

#' Alias for descriptivesByCondition
#' @export
TSDescriptivesByCondition <- descriptivesByCondition

#' Plot Time Series Descriptive Statistics
#'
#' @param data Data frame of descriptive statistics.
#' @param title Plot title.
#' @param vars Optional variables to melt.
#' @return A ggplot object.
#' @export
plotTSDescriptives <- function(data,title="",vars=NULL) {
  if(is.null(vars)){
    data <- melt(as.data.frame(data),id.vars = c("Start","End","Condition"), )
  }
  else {
    data <- melt(as.data.frame(data),id.vars = c("Start","End","Condition"),measure.vars = vars )
    
  }
  #levels(data$Condition) <- unique(data$Condition[order(data$Start)])
  data$value[!is.finite(data$value)] <- NA
  
  y.max <- max(data$value, na.rm=T)
  y.min <- min(data$value, na.rm=T)
  yrange <- y.max - y.min
  rect.max <- (y.min - (yrange/20.0))
  rect.min <- rect.max - (yrange/20.0)
  g1<- ggplot(data = data,mapping = aes(x=(Start + (End-Start)/2),y=value,col=variable)) + geom_line(size=1) + geom_point(size=3) + xlab("Time") + ylab("Parameter Value") + ggtitle(title)
  g1 <- g1 + geom_rect(aes_string(xmin="Start",xmax="End",ymin=rect.min,ymax=rect.max, fill="Condition",col=NULL))
  g1 <- g1 + guides(colour = guide_legend(title="Parameter",ncol = 1,override.aes = list(fill=NA)), fill = guide_legend(ncol = 1)) + theme(legend.position = "right")

  g1 <- g1 + scale_y_continuous()
  g1 <- g1 + scale_x_datetime()

  return(g1)
}

#' Interrupted Time Series Analysis (Not Implemented)
#'
#' @param data Data frame or time series data.
#' @param col.data Column name for the data (default "EDA").
#' @param col.condition Column name for the condition (default "Condition").
#' @return Not implemented.
#' @export
interruptedTimeSeries <- function(data, col.data="EDA", col.condition="Condition") {
}

#' Calculate p-value from an lm object
#'
#' @param modelobject An object of class 'lm' or 'systemfit.equation'.
#' @return The p-value of the F-statistic.
#' @export
lmp <- function (modelobject) {
  if (!(class(modelobject) == "lm" || class(modelobject) == "systemfit.equation")) stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#' Calculate Euclidean Distance
#'
#' @param x Numeric value.
#' @param y Numeric value.
#' @return Euclidean distance computed as sqrt(x*x + y*y).
#' @export
distance <- function(x,y){
  return(sqrt(x*x + y*y));
}

#' Visualize Vector Field
#'
#' @param x Numeric vector for x-axis.
#' @param y Numeric vector for y-axis.
#' @param f Function for vector field computation.
#' @param scale Scaling factor (default 0.1).
#' @export
vectorflow <- function(x,y,f,scale=0.1) {
  xx <- c(min(x),max(x))
  yy <- c(min(y),max(y))
  vectorfield(f,xx,yy,scale=scale)

}

#' Plot Topography Density
#'
#' @param x Numeric vector for x-axis.
#' @param y Numeric vector for y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @return A ggplot object showing 2D density plot.
#' @export
topography <- function(x,y,xlab="X",ylab="Y") {
  subdata <- as.data.frame(x)
  subdata$x <- x
  subdata$y <- y
  ggplot(subdata,aes(x=x,y=y))+geom_density2d()+xlab(xlab)+ylab(ylab)+xlim(min(x),max(x))+ylim(min(y),max(y))
}

#' Normalize a Numeric Vector
#'
#' @param x Numeric vector.
#' @return Normalized vector where the mean is subtracted and divided by the standard deviation.
#' @export
normalize <- function(x){
  return((x-mean(x,na.rm = T))/sd(x,na.rm = T))
}

#' Create a Vector Flow Plot
#'
#' @param x Numeric vector for x coordinates.
#' @param y Numeric vector for y coordinates.
#' @param xlabel Label for x-axis.
#' @param ylabel Label for y-axis.
#' @param title Plot title.
#' @return A ggplot object displaying the vector flow plot.
#' @export
vectorfield <- function(x,y,xlabel="X",ylabel="Y",title="Vector Flow Plot"){
  subdata <- as.data.frame(x)
  subdata$x <- x
  subdata$y <- y
  subdata$x.lag1 <- Lag(x,1)
  subdata$x.lead1 <- Lag(x,-1)
  subdata$x.vel <- (subdata$x.lead1 - subdata$x.lag1)/2
  subdata$y.lag1 <- Lag(y,1)
  subdata$y.lead1 <- Lag(y,-1)
  subdata$y.vel <- (subdata$y.lead1 - subdata$y.lag1)/2
  #Graph of a Vector Plot
  ggplot(data=subdata, aes(x=x,y=y))+geom_density2d()+geom_segment(aes(xend=x+x.vel, yend=y+y.vel),arrow=arrow(length=unit(.2,"cm")))+labs(list(title=title, x= xlabel, y = ylabel))+xlim(min(x,na.rm = T),max(x,na.rm = T))+ylim(min(y,na.rm = T),max(y,na.rm = T))

}

#' Simulate Data using Linear Models
#'
#' @param xmodel Linear model for x simulation.
#' @param ymodel Linear model for y simulation.
#' @param n Number of samples to simulate (default 1000).
#' @param start Starting values as a vector (default c(1,2)).
#' @return A matrix with simulated data (n rows x 2 columns).
#' @export
simulate <- function(xmodel,ymodel,n=1000, start=c(1,2)) {
  x = start[1]
  y = start[2]
  out = matrix(data=NA,nrow = n, ncol = 2)
  out[1,1] <- x
  out[1,2] <- y
  for(i in seq(2,n)){
    dx <- predict(xmodel,data.frame(x=out[i-1,1],y=out[i-1,2]))[[1]]
    dy <- predict(xmodel,data.frame(x=out[i-1,1],y=out[i-1,2]))[[1]]
    x = out[i-1,1] + dx
    y = out[i-1,2] + dy
    out[i,1] = x
    out[i,2] = y
  }
  return(out)

}

#' Build Dyadic Model
#'
#' @param p1 Data for participant 1.
#' @param p2 Data for participant 2.
#' @param type Model type (default 2).
#' @param step Step size (default 0.25).
#' @param measure Column name for measurement (default "EDA").
#' @return A model object.
#' @export
dyad.dsmodel <- function(p1,p2,type=2, step=0.25, measure="EDA") {
  p1.name <- deparse(substitute(p1))
  p2.name <- deparse(substitute(p2))
  data <- as.dyad(p1,p2, cols=c(measure))

  return(dsmodel(data[,2],data[,3],type=type, step=step))
}

#' Plot Dyadic Graph
#'
#' @param p1 Data for participant 1.
#' @param p2 Data for participant 2.
#' @param norm Logical indicating whether to normalize the data (default TRUE).
#' @param epoch Epoch size (default 1).
#' @param step Step size (default 0.25).
#' @param title Plot title.
#' @param type Model type (default 2).
#' @param measure Column name for measurement (default "EDA").
#' @return A plot of the model params
#' @export
dyad.dsgraph <- function(p1,p2,norm=T,epoch=1, step=0.25, title="State Space Plot", type=2,measure="EDA" ) {
  p1.name <- deparse(substitute(p1))
  p2.name <- deparse(substitute(p2))
  data <- as.dyad(p1,p2,cols=c(measure))
  x <- data[,2]
  y <- data[,3]
  downsample_factor <- as.integer(round(epoch/as.numeric(mean(diff(data$Timestamp)))))
  print(paste("Downsampling both signals by",downsample_factor))
  x<- decimate(x,q = downsample_factor)
  y<- decimate(y,q = downsample_factor)

  statespacegraph(x,y,norm=norm, step=step,title=title, xlabel=p1.name, ylabel=p2.name,type=2)

}

#' Compute Lag Difference
#'
#' @param x Numeric vector.
#' @param lag Lag value.
#' @param robust Logical indicating use of robust method (default TRUE).
#' @return Lagged difference.
#' @export
o.Lag <- function(x,lag=0, robust=T) {
  if(robust){
    dx <- (Lag(x,(-1*lag)-1) - Lag(x,(-1*lag)+1))/2
  }
  else {
    dx <- Lag(x,(-1*lag)-1) -  Lag(x,(-1*lag))
  }
  
  return(dx)
}

#' Estimate Lag with Moving Window
#'
#' @param x Numeric vector.
#' @param window Window size (default 10).
#' @return Estimated lag values.
#' @export
o.estLag <- function(x,window=10){
  mat <- matrix(nrow=length(x),ncol = window)
  for(n in 1:window){
    l <- n - window/2
    mat[,n] <- o.Lag(x,l)
  }
  return(rowWeightedMeans(mat,w=gausswin(window),na.rm = T))
}

#' Smooth a Numeric Vector
#'
#' @param x Numeric vector.
#' @param window Window size (default 10).
#' @return Smoothed vector.
#' @export
o.smooth <- function(x,window=10){
  mat <- matrix(nrow=length(x),ncol = window)
  for(n in 1:window){
    l <- n - window/2
    mat[,n] <- Lag(x,l)
  }
  return(rowWeightedMeans(mat,w=gausswin(window),na.rm = T))
}

#' Generate Random Dyadic Data
#'
#' @param dur Duration in seconds (default 60*120).
#' @param fs Sampling rate (default 32).
#' @return Data frame with Timestamp and two fake signals.
#' @export
rand.dyad <- function(dur=60*120,fs=32){
  now = Sys.time()
  l = dur*fs
  ts <- seq(now,(now+as.difftime(dur,unit="secs")),as.difftime(1.0/fs,units="secs"))[1:l]
  fake.1 <- (sample(0:10000,l,replace = T)/1000.0)[1:l]
  fake.2 <- (sample(0:10000,l,replace = T)/1000.0)[1:l]
  d <- data.frame(Timestamp=ts,fake.1=fake.1,fake.2=fake.2)
  return(d)
}

#' Extract Coefficient with Significance Check
#'
#' @param mdl Linear model object.
#' @param c Parameter index (default 2).
#' @param cutoff Significance cutoff (default 0.05).
#' @param useD Logical, if TRUE use standardized coefficient (default FALSE).
#' @return Coefficient value or 0 if not significant.
#' @export
o.coef <- function(mdl,c=2,cutoff=0.05,useD=F){
  s <- summary(mdl)
  p <- s$coefficients[c,4]
  if(p > cutoff){
    return(0)
  }
  else if(useD) {
    return(s$coefficients[c,1]/s$coefficients[c,2])
  }
  else {
    return(s$coefficients[c,1])
  }
}

#' Fit Dyadic Model
#'
#' @param x Numeric vector for signal x.
#' @param y Numeric vector for signal y.
#' @param type Model type (default 2).
#' @param step Step size (default 0.25).
#' @param p.value P-value threshold (default 0.01).
#' @param verbose Logical indicating verbose output (default FALSE).
#' @return A list containing the fitted models, R-squared values, prediction functions, and model equations.
#' @export
dsmodel <- function(x,y,type=2, step=0.25,p.value=0.01,verbose=F){
  x_prime <-o.Lag(x)
  y_prime <- o.Lag(y)
  base_x_model <- lm(x_prime ~ x)
  base_y_model <- lm(y_prime ~ y)
  
  s1 <- ""
  s2 <- ""
  
  if(type==1){
    x_model <- lm(x_prime ~ x * y)
    y_model <- lm(y_prime ~ y * x)

    b0 <- coef(x_model)[[1]]
    b1 <- coef(x_model)[[2]]
    b2 <- coef(x_model)[[3]]
    bxy <- coef(x_model)[[4]]
    b3 <- coef(y_model)[[1]]
    b4 <- coef(y_model)[[2]]
    b5 <- coef(y_model)[[3]]
    s1 <- cat("x' = ",b0," + ",b1,"*x"," + ",b2,"*y","\n")
    s2 <- cat("y' = ",b3," + ",b4,"*y"," + ",b5,"*x","\n")
  }
  else if(type == 2){
    x_model <- lm(x_prime ~ x + I(y-x))
    y_model <- lm(y_prime ~ y + I(x-y))

    b0 <- o.coef(x_model,1)
    b1 <- o.coef(x_model,2)
    b2 <- o.coef(x_model,3)
    b3 <- o.coef(y_model,1)
    b4 <- o.coef(y_model,2)
    b5 <- o.coef(y_model,3)
    s1 <- paste("x'"," = ",b0," + ",b1,"x"," + ",b2,"(y-x)",sep="")
    s2 <- paste("y'"," = ",b3," + ",b4,"*y"," + ",b5,"(x-y)",sep="")
  }
  else if (type == 4) {
      x_model <- lm(x_prime ~ I(mean(x,na.rm=T)-x) * I(y-x))
      y_model <- lm(y_prime ~ I(mean(y,na.rm=T)-y) * I(x-y) )
      b0 <- o.coef(x_model,1)
      b1 <- o.coef(x_model,2)
      b2 <- o.coef(x_model,3)
      b3 <- o.coef(y_model,1)
      b4 <- o.coef(y_model,2)
      b5 <- o.coef(y_model,3)
      s1 <- cat("x' = ",b0," + ",b2,"*(y-x)","\n")
      s2 <- cat("y' = ",b3," + ",b5,"*(x-y)","\n")
  }
  else if (type == 0) {
    x_model <- lm(x_prime ~ x)
    y_model <- lm(y_prime ~ y)
    b0 <- coef(x_model)[[1]]
    b1 <- 0
    b2 <- coef(x_model)[[2]]
    b3 <- coef(y_model)[[1]]
    b4 <- 0
    b5 <- coef(y_model)[[2]]
    s1 <- cat("x' = ",b0," + ",b2,"*(y-x)","\n")
    s2 <- cat("y' = ",b3," + ",b5,"*(x-y)","\n")
  }
else if (type == 5) {
  x_model <- lm(x_prime ~  I(y-x))
  y_model <- lm(y_prime ~ I(x-y) )
  b0 <- o.coef(x_model,1)
  b1 <- 0
  b2 <- o.coef(x_model,2)
  b3 <- o.coef(y_model,1)
  b4 <- 0
  b5 <- o.coef(y_model,2)
  s1 <- cat("x' = ",b0," + ",b2,"*(y-x)","\n")
  s2 <- cat("y' = ",b3," + ",b5,"*(x-y)","\n")
}
  else {
    x_model <- lm(x_prime ~ x + y)
    y_model <- lm(y_prime ~ y + x)

    b0 <- coef(x_model)[[1]]
    b1 <- coef(x_model)[[2]]
    b2 <- coef(x_model)[[3]]
    b3 <- coef(y_model)[[1]]
    b4 <- coef(y_model)[[2]]
    b5 <- coef(y_model)[[3]]
    s1 <- cat("x' = ",b0," + ",b1,"*x"," + ",b2,"*y","\n")
    s2 <- cat("y' = ",b3," + ",b4,"*y"," + ",b5,"*x","\n")
  }
  x_f <- function(xt,yt){
    val <- 0
    try(val <- predict(object = x_model,newdata=data.frame(x=xt,y=yt))[[1]])
    return(val)
  }
  y_f <- function(xt,yt){
    val <- 0
    try(val <- predict(object = y_model,newdata=data.frame(x=xt,y=yt))[[1]])
    return(val)
  }
  xx <- c(min(x),max(x))
  yy <- c(min(y),max(y))

  if(verbose){
    print(s1)
    print(summary(x_model))
    print(s2)
    print(summary(y_model))
    
  }

  xvals <- seq(min(x,na.rm = T),max(x,na.rm = T),step)
  yvals <- seq(min(y,na.rm = T),max(y,na.rm = T),step)
  subdata <- matrix(nrow = length(xvals)*length(yvals), ncol=4)
  n= 0
  for(xt in xvals){
    for(yt in yvals){
      subdata[n,] <- c(xt,yt,x_f(xt,yt),y_f(xt,yt))
      n= n + 1
    }
  }

  colnames(subdata) <- c("x","y","xprime","yprime")
  subdata <- as.data.frame(subdata)
    dx.r.squared <- 0
    x.r.squared <- 0
    x.base.r.squared <- 0
    
    if(lmp(x_model) < p.value){
      x.r.squared <- summary(x_model)$adj.r.squared
      x.base.r.squared <- summary(base_x_model)$adj.r.squared
      dx.r.squared <- x.r.squared - x.base.r.squared
    }
  
    y.r.squared <- 0
    y.base.r.squared <- 0
    dy.r.squared <- 0
  
    if(lmp(y_model) < p.value){
      y.r.squared <- summary(y_model)$adj.r.squared
      y.base.r.squared <- summary(base_y_model)$adj.r.squared
      dy.r.squared <- y.r.squared - y.base.r.squared
    }
  
  return(list(x_model=x_model,y_model=y_model,x.base.r.squared=x.base.r.squared,y.base.r.squared=y.base.r.squared, x.r.squared=x.r.squared, y.r.squared=y.r.squared,dx.r.squared=dx.r.squared,dy.r.squared=dy.r.squared, mat=subdata,x_func=x_f, y_func=y_f,x.eq=s1,y.eq=s2,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,b5=b5))
}

#' Fit Dyadic Model using SEM
#'
#' @param d Data frame containing dyad data.
#' @param type Model type (default 2).
#' @param step Step size (default 0.25).
#' @param p.value P-value threshold (default 0.01).
#' @param verbose Logical indicating verbose output (default TRUE).
#' @return A SEM model object.
#' @export
dsmodel.sem <- function(d,type=2, step=0.25,p.value=0.01,verbose=T){
  x <- d[,2]
  y <- d[,3]
  x_prime <-o.Lag(x)
  y_prime <- o.Lag(y)
  dxy <- (x-y)
  dyx <- (y-x)
  
  s1 <- ""
  s2 <- ""
  
  data <- data.frame(x_prime=x_prime,y_prime=y_prime,x=x,y=y,dxy=dxy,dyx=dyx)
    model <- '
  # Regression model. The b1 and b2 are labels
  # similar to (b1) and (b2) in Mplus. These are
  # needed to perform wald test.
  # jobperf ~ b1*wbeing + b2*jobsat
    Reg.X =~ a1*x + a2*dyx
    Reg.Y =~ b1*y + b2*dyx
    x_prime ~ c1*Reg.X + c2*Reg.Y
    y_prime ~ d1*Reg.Y + d2*Reg.X
  # Variances
    x ~~ x
    y ~~ y
    dyx ~~ dyx
    dxy ~~ dxy
  
  # Covariance/correlation
    x ~~ y
  '
    
    mdl <- sem(model,data=data,missing = "fiml")
    return(mdl)
}

#' Fit Multi-lag Dyadic Model
#'
#' @param x Numeric vector for signal x.
#' @param y Numeric vector for signal y.
#' @param step Step size (default 0.25).
#' @param p.value P-value threshold (default 0.01).
#' @param type Model type (default 2).
#' @return A list containing models, statistics, and prediction functions.
#' @export
dsmodel.multilag <- function(x,y, step=0.25,p.value=0.01,type=2){
  x_prime <-o.Lag(x)
  y_prime <- o.Lag(y)
  s1 <- ""
  s2 <- ""
  base_x_model <- lm(x_prime ~ x)
  base_y_model <- lm(y_prime ~ y)
    x.lag1 <- Lag(x, -1)
    y.lag1 <- Lag(y,-1)
    x.lag2 <- Lag(x, -2)
    y.lag2 <- Lag(y,-2)
    x.lag3 <- Lag(x, -3)
    y.lag3 <- Lag(y,-3)
  
    x_model <- lm(x_prime ~ I(y-x) + I(y.lag1-x.lag1) + I(y.lag2-x.lag2) + I(y.lag3-x.lag3))
    y_model <- lm(y_prime ~ I(x-y) + I(x.lag1 - y.lag1) + I(x.lag2 - y.lag2)+I(x.lag3 - y.lag3))
    
    b0 <- coef(x_model)[[1]]
    b1 <- coef(x_model)[[2]]
    b2 <- coef(x_model)[[3]]
    b3 <- coef(y_model)[[1]]
    b4 <- coef(y_model)[[2]]
    b5 <- coef(y_model)[[3]]
    s1 <- paste("x'"," = ",b0," + ",b1,"x"," + ",b2,"(y-x)",sep="")
    s2 <- paste("y'"," = ",b3," + ",b4,"*y"," + ",b5,"(x-y)",sep="")

  x_f <- function(xt,yt){
    val <- predict(object = x_model,newdata=data.frame(x=xt,y=yt))[[1]]
    return(val)
  }
  y_f <- function(xt,yt){
    val <- predict(object = y_model,newdata=data.frame(x=xt,y=yt))[[1]]
    return(val)
  }
  xx <- c(min(x),max(x))
  yy <- c(min(y),max(y))
  
  print(s1)
  print(summary(x_model))
  print(s2)
  print(summary(y_model))
  
  xvals <- seq(min(x,na.rm = T),max(x,na.rm = T),step)
  yvals <- seq(min(y,na.rm = T),max(y,na.rm = T),step)
  subdata <- matrix(nrow = length(xvals)*length(yvals), ncol=4)
  n= 0
  # Note: loop for filling subdata is commented out in original code
  subdata <- as.data.frame(subdata)
  dx.r.squared <- 0
  x.r.squared <- 0
  x.base.r.squared <- 0
  
  if(lmp(x_model) < p.value){
    x.r.squared <- summary(x_model)$r.squared
    x.base.r.squared <- summary(base_x_model)$r.squared
    dx.r.squared <- (x.r.squared - x.base.r.squared)
  }
  
  y.r.squared <- 0
  y.base.r.squared <- 0
  dy.r.squared <- 0
  
  if(lmp(y_model) < p.value){
    y.r.squared <- summary(y_model)$r.squared
    y.base.r.squared <- summary(base_y_model)$r.squared
    dy.r.squared <- (y.r.squared - y.base.r.squared)
  }
  
  return(list(x_model=x_model,y_model=y_model,x.base.r.squared=x.base.r.squared,y.base.r.squared=y.base.r.squared, x.r.squared=x.r.squared, y.r.squared=y.r.squared,dx.r.squared=dx.r.squared,dy.r.squared=dy.r.squared, mat=subdata,x_func=x_f, y_func=y_f,x.eq=s1,y.eq=s2,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,b5=b5))
}

#' Rescale a Vector to [0,1]
#'
#' @param x Numeric vector.
#' @return The rescaled vector.
#' @export
o.rescale <- function(x){
  out <- (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))
  return(out)
}

#' Plot State Space Graph
#'
#' @param xt Numeric vector for x values.
#' @param yt Numeric vector for y values.
#' @param norm Logical to normalize data (default TRUE).
#' @param step Step size (default 0.25).
#' @param title Plot title.
#' @param xlabel Label for x-axis.
#' @param ylabel Label for y-axis.
#' @param type Model type (default 3).
#' @return A ggplot object representing the state space graph.
#' @export
statespacegraph <- function(xt,yt,norm=T, step=0.25,title="State Space Plot", xlabel="X Data", ylabel="Y Data",type=3) {
  xlab <- xlabel
  ylab <- ylabel
  if(norm){
    xt <- o.rescale(xt)
    yt <- o.rescale(yt)
    xlabel <- paste(xlabel,"(Std.)")
    ylabel <- paste(ylabel,"(Std.)")
  }
  s <- dsmodel(xt,yt,step=step,type=type)
  data <- s$mat

  p <- ggplot(data=data, aes(x=x,y=y))+geom_segment(aes(xend=(x+xprime), yend=(y+yprime),lineend="round", colour=(sqrt(xprime*xprime+yprime*yprime))),arrow=arrow(length=unit(.2,"cm")))+
    scale_colour_gradient(low="#0000CC",high="#CC0000",name="Distance")+
    labs(list(title=title, x=xlabel, y =ylabel))+
    xlim(min(xt,na.rm = T),max(xt,na.rm = T))+
    ylim(min(yt,na.rm = T),max(yt,na.rm = T))
  lab.1 <- expression(paste(R^2," ",xlabel,"="))
  lab.2 <- expression(paste(R^2," ",ylabel,"="))
  r2.x <- round(summary(s$x_model)$r.squared,digits = 2)
  r2.y <- round(summary(s$y_model)$r.squared,digits = 2)
  b2 <- round(standardCoefs(s$x_model)[2,2],digits=2)
  b5 <- round(standardCoefs(s$y_model)[2,2],digits=2)
  footnote <- substitute(paste(R^2," ",xlabel,"=", r2.x," | ",R^2," ",ylabel,"=", r2.y, " | ",beta[xlab],"=",b2," | ",beta[ylab],"=",b5))
  g <- arrangeGrob(p, sub = textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 16)))

  g
}

#' Create Plot for Dyadic State Space Graph
#'
#' @param s A dyadic statespace model object.
#' @param title Plot title.
#' @param xlabel Label for x-axis.
#' @param ylabel Label for y-axis.
#' @return A ggplot object.
#' @export
plot.dsgraph <- function(s,title="State Space Plot", xlabel="X Data", ylabel="Y Data") {
  xlab <- xlabel
  ylab <- ylabel
  data <- s$mat
  xt <- s$mat$x
  yt <- s$mat$y
  p <- ggplot(data=data, aes(x=x,y=y))+geom_segment(aes(xend=(x+xprime), yend=(y+yprime),lineend="round", colour=(sqrt(xprime*xprime+yprime*yprime))),arrow=arrow(length=unit(.2,"cm")))+
    scale_colour_gradient(low="#0000CC",high="#CC0000",name="Distance")+
    labs(list(title=title, x=xlabel, y =ylabel))+
    xlim(min(xt,na.rm = T),max(xt,na.rm = T))+
    ylim(min(yt,na.rm = T),max(yt,na.rm = T))
  return(p)
}

#' Apply Moving Window Operation
#'
#' @param x Numeric vector.
#' @param window_size Size of each window.
#' @param window_overlap Overlap between windows.
#' @param FUN Function to apply on each window.
#' @param na.rm Logical indicating whether to remove NAs (default TRUE).
#' @return A vector of computed values.
#' @export
window <- function(x, window_size, window_overlap=0, FUN, na.rm=T) {
  x <- as.vector(x)
  output_length = length(x)/(window_size- window_overlap/2)
  output = vector(length=output_length)
  n = 0
  step = window_size-window_overlap
  for(i in seq(1,length(x), step)){
    start = i
    end = i + step
    if(end > length(x)){
      end = length(x)
    }
    output[n] <- FUN(x[start:end])
    n = n+1
  }
  
  return(output)
}

#' Copy Text to Clipboard
#'
#' @param x Object to be copied.
#' @export
pbcopy <- function(x){
  write.table(file = pipe("pbcopy"), x, sep = "\t",quote=F)
}
