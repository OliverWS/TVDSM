  
require(signal)
require(matrixStats)
require(Hmisc)
library(lavaan)

RMSSD <- function(hr.period) {
  hr.period.delta <- diff(hr.period)
  return(sqrt(mean(hr.period.delta*hr.period.delta)))
}

linearTrend <- function(data) {
  t <- seq_along(data)
  mdl <- lm(data ~ t)
  trend <- coef(mdl)[[2]]
  return(trend)
}
se <- function(x) sqrt(var(x)/length(x))

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


dummyCodeByCondition <- function(filename,codes,outputFile=paste(filename,"_DummyCoded.csv",sep=""),saveFile=F) {
  
  if(typeof(codes) == "character"){
    conditions <- read.Codes(codes)
  }
  else {
    conditions <- codes
  }
  eda <- read.eda(file = filename)
  fs <- getFS(eda)
  nConditions <- dim(conditions)[1]
  tz <- attr(eda$Timestamp[1],"tz")
  eda$Condition <- NA
  for (i in 1:nConditions) {
    start <- conditions[["Start Time"]][i]
    end <- conditions[["End Time"]][i]
    condition <- as.character(conditions[["Condition"]][i])
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
    write.csv(eda,file = outputFile,quote = T)

  }
  return(eda)
}


descriptivesByCondition <- function(filename, codes,col="EDA",title=filename) {
  
  if(typeof(codes) == "character"){
    conditions <- read.Codes(codes)
  }
  else {
    conditions <- codes
  }
  eda <- read.eda(file = filename)
  fs <- getFS(eda)
  nConditions <- dim(conditions)[1]
  results <- list()
  tz <- attr(eda$Timestamp[1],"tz")
  for (i in 1:nConditions) {
    start <- conditions[["Start Time"]][i]
    end <- conditions[["End Time"]][i]
    condition <- conditions[["Condition"]][i]
    validTimes <- ((eda$Timestamp >= start) & (eda$Timestamp < end))
    data <- subset(eda,validTimes)
    desc <- tsDescriptives(data[[col]])
    
    results[[i]] <- desc
  }

  results <- as.data.frame(t(as.data.frame(sapply(results, unlist))))
  rownames(results) <- NULL
  #Fix weird issue where values get turned into text
  results <- lapply(results,as.double)
  results$Condition <- conditions$Condition
  results$Start <- conditions$`Start Time`
  results$End <- conditions$`End Time`
  
  #g2 <-ggplot(data = subset(eda, ((eda$Timestamp >= min(results$Start)) & (eda$Timestamp < max(results$End)))) ,mapping = aes(x=Timestamp,y=EDA)) + geom_line(size=1,col="#1FBFC4") + xlab("Time") + ylab("EDA")

  g1 <- plotTSDescriptives(results,title=title)
  
  #plt <- plot_grid(g1,g2,ncol=1,nrow=2,align = "h")
  print(g1)
  return(results)
  
  
}

TSDescriptivesByCondition <- descriptivesByCondition
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


interruptedTimeSeries <- function(data, col.data="EDA", col.condition="Condition") {
}




lmp <- function (modelobject) {
  if (!(class(modelobject) == "lm" || class(modelobject) == "systemfit.equation")) stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}






distance <- function(x,y){
  return(sqrt(x*x + y*y));
}

vectorflow <- function(x,y,f,scale=0.1) {
  xx <- c(min(x),max(x))
  yy <- c(min(y),max(y))
  vectorfield(f,xx,yy,scale=scale)

}

topography <- function(x,y,xlab="X",ylab="Y") {
  subdata <- as.data.frame(x)
  subdata$x <- x
  subdata$y <- y
  ggplot(subdata,aes(x=x,y=y))+geom_density2d()+xlab(xlab)+ylab(ylab)+xlim(min(x),max(x))+ylim(min(y),max(y))
}

normalize <- function(x){
  return((x-mean(x,na.rm = T))/sd(x,na.rm = T))
}
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
dyad.statespace <- function(p1,p2,type=2, step=0.25, measure="EDA") {
  p1.name <- deparse(substitute(p1))
  p2.name <- deparse(substitute(p2))
  data <- as.dyad(p1,p2, cols=c(measure))

  return(statespace(data[,2],data[,3],type=type, step=step))
}

dyad.statespacegraph <- function(p1,p2,norm=T,epoch=1, step=0.25, title="State Space Plot", type=2,measure="EDA" ) {
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

o.Lag <- function(x,lag=0, robust=F) {
  if(robust){
    dx <- (Lag(x,(-1*lag)-1) - Lag(x,(-1*lag)+1))/2
  }
  else {
    dx <- Lag(x,(-1*lag)-1) -  Lag(x,(-1*lag))
  }
  
  return(dx)
}


o.estLag <- function(x,window=10){
  mat <- matrix(nrow=length(x),ncol = window)
  for(n in 1:window){
    l <- n - window/2
    mat[,n] <- o.Lag(x,l)
  }
  return(rowWeightedMeans(mat,w=gausswin(window),na.rm = T))
  
}

o.smooth <- function(x,window=10){
  mat <- matrix(nrow=length(x),ncol = window)
  for(n in 1:window){
    l <- n - window/2
    mat[,n] <- Lag(x,l)
  }
  return(rowWeightedMeans(mat,w=gausswin(window),na.rm = T))
  
}
rand.dyad <- function(dur=60*120,fs=32){
  now = Sys.time()
  l = dur*fs
  ts <- seq(now,(now+as.difftime(dur,unit="secs")),as.difftime(1.0/fs,units="secs"))[1:l]
  fake.1 <- (sample(0:10000,l,replace = T)/1000.0)[1:l]
  fake.2 <- (sample(0:10000,l,replace = T)/1000.0)[1:l]
  d <- data.frame(Timestamp=ts,fake.1=fake.1,fake.2=fake.2)
  return(d)
}

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

statespace <- function(x,y,type=2, step=0.25,p.value=0.01,verbose=F){
  #(Xt+1 – Xt)= b0+b1(Xt)+b2(Yt)+e
  #(Yt+1 – Yt)= b3+b4(Yt)+b5(Xt)+e
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
    byx <- coef(y_model)[[4]]
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
#     x_prime <- o.estLag(x,window = 32)
#     y_prime <- o.estLag(y,window=32)
#     base_x_model <- lm(x_prime ~ x)
#     base_y_model <- lm(y_prime ~ y)
#     
#     dyx <- o.smooth(y-x,window = 32)
#     dxy <- o.smooth(x-y,window = 32)
#     
    
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
  #     x_prime <- o.estLag(x,window = 32)
  #     y_prime <- o.estLag(y,window=32)
  #     base_x_model <- lm(x_prime ~ x)
  #     base_y_model <- lm(y_prime ~ y)
  #     
  #     dyx <- o.smooth(y-x,window = 32)
  #     dxy <- o.smooth(x-y,window = 32)
  #     
  
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



library(systemfit)

statespace.sem <- function(d,type=2, step=0.25,p.value=0.01,verbose=T){
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
#     x_model <- fit$eq[[1]]
#     y_model <- fit$eq[[2]]
#     
#     b0 <- o.coef(x_model,1)
#     b1 <- o.coef(x_model,2)
#     b2 <- o.coef(x_model,3)
#     b3 <- o.coef(y_model,1)
#     b4 <- o.coef(y_model,2)
#     b5 <- o.coef(y_model,3)
# 
# 
#   if(verbose){
#     
#     print(summary(x_model))
#     print(summary(y_model))
#     
#     
#   }
#   
#   
#   x.base.r.squared=0
#   y.base.r.squared=0
#   x.r.squared = summary(x_model)$adj.r.squared
#   y.r.squared = summary(y_model)$adj.r.squared
#   
  
  #return(list(x_model=x_model,y_model=y_model,model=fit,x.base.r.squared=x.base.r.squared,y.base.r.squared=y.base.r.squared, x.r.squared=x.r.squared, y.r.squared=y.r.squared,dx.r.squared=x.r.squared,dy.r.squared=y.r.squared,x.eq=s1,y.eq=s2,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,b5=b5))
  return(mdl)
}





statespace.multilag <- function(x,y, step=0.25,p.value=0.01,type=2){
  
  #(Xt+1 – Xt)= b0+b1(Xt)+b2(Yt)+e
  #(Yt+1 – Yt)= b3+b4(Yt)+b5(Xt)+e
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
#   for(xt in xvals){
#     for(yt in yvals){
#       subdata[n,] <- c(xt,yt,x_f(xt,yt),y_f(xt,yt))
#       n= n + 1
#     }
#   }
#   
#   colnames(subdata) <- c("x","y","xprime","yprime")
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





o.rescale <- function(x){
  out <- (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))
  return(out)
}

statespacegraph <- function(xt,yt,norm=T, step=0.25,title="State Space Plot", xlabel="X Data", ylabel="Y Data",type=3) {
  xlab <- xlabel
  ylab <- ylabel
  if(norm){
    xt <- o.rescale(xt)
    yt <- o.rescale(yt)
    xlabel <- paste(xlabel,"(Std.)")
    ylabel <- paste(ylabel,"(Std.)")

  }
  s <- statespace(xt,yt,step=step,type=type)
  data <- s$mat

  p <- ggplot(data=data, aes(x=x,y=y))+geom_segment(aes(xend=(x+xprime), yend=(y+yprime),lineend="round", colour=(sqrt(xprime*xprime+yprime*yprime))),arrow=arrow(length=unit(.2,"cm")))+scale_colour_gradient(low="#0000CC",high="#CC0000",name="Distance")+labs(list(title=title, x=xlabel, y =ylabel))+xlim(min(xt,na.rm = T),max(xt,na.rm = T))+ylim(min(yt,na.rm = T),max(yt,na.rm = T))
  # source: http://bigdata-analyst.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
  lab.1 <- expression(paste(R^2," ",xlabel,"="))
  lab.2 <- expression(paste(R^2," ",ylabel,"="))
  r2.x <- round(summary(s$x_model)$r.squared,digits = 2)
  r2.y <- round(summary(s$y_model)$r.squared,digits = 2)
  b2 <- round(standardCoefs(s$x_model)[2,2],digits=2)
  b5 <- round(standardCoefs(s$y_model)[2,2],digits=2)
  footnote <- substitute(paste(R^2," ",xlabel,"=", r2.x," | ",R^2," ",ylabel,"=",r2.y, " | ",beta[xlab],"=",b2," | ",beta[ylab],"=",b5))
  g <- arrangeGrob(p, sub = textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 16)))

  g

}

plot.statespace <- function(s,title="State Space Plot", xlabel="X Data", ylabel="Y Data") {
  xlab <- xlabel
  ylab <- ylabel
  data <- s$mat
  xt <- s$mat$x
  yt <- s$mat$y
  p <- ggplot(data=data, aes(x=x,y=y))+geom_segment(aes(xend=(x+xprime), yend=(y+yprime),lineend="round", colour=(sqrt(xprime*xprime+yprime*yprime))),arrow=arrow(length=unit(.2,"cm")))+scale_colour_gradient(low="#0000CC",high="#CC0000",name="Distance")+labs(list(title=title, x=xlabel, y =ylabel))+xlim(min(xt,na.rm = T),max(xt,na.rm = T))+ylim(min(yt,na.rm = T),max(yt,na.rm = T))
  # source: http://bigdata-analyst.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
#   lab.1 <- expression(paste(R^2," ",xlabel,"="))
#   lab.2 <- expression(paste(R^2," ",ylabel,"="))
#   r2.x <- round(summary(s$x_model)$r.squared,digits = 2)
#   r2.y <- round(summary(s$y_model)$r.squared,digits = 2)
#   b2 <- round(standardCoefs(s$x_model)[2,2],digits=2)
#   b5 <- round(standardCoefs(s$y_model)[2,2],digits=2)
#   footnote <- substitute(paste(R^2," ",xlabel,"=", r2.x," | ",R^2," ",ylabel,"=",r2.y, " | ",beta[xlab],"=",b2," | ",beta[ylab],"=",b5))
  return(p)
  
}



#
# lag.matrix <- function(x1, x2, func=lag.range=c(-10,10)) {
#   nlags = max(lag.range) - min(lag.range)
#   mat = matrix(data=0, ncol = length(x1), nrow = nlags)
#   data <- cbind(x1,x2)
#   for(lag in seq(from = min(lag.range), to = max(lag.range),by = 1)){
#
#   }
# }


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





pbcopy <- function(x){
  write.table(file = pipe("pbcopy"), x, sep = "\t",quote=F)
}
