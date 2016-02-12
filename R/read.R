library(readxl)
read.ERCodes <- function(path){
  if(grepl(".xls", path)){
    data <- read_excel(path,"text")
  }
  else {
    data <- read.csv(path,header = T)
  }
  data["Time.Start"] <- as.difftime(as.character(data$Time.Start),format = "%M:%S",units = "secs")
  data["Time.End"] <- as.difftime(as.character(data$Time.End),format = "%M:%S",units = "secs")
  output <- data.frame(Timestamps=(0:as.numeric(max(data$Time.End),units="secs")))
  for(t in output$Timestamps)
  {
    t <- which(data$Time.Start <= t)
    idx <- which((data$Time.Start >= t) & (data$Time.End <= t))
    state <- data$Engage.State[idx]
    if(!is.null(state)){
      output$Engage.State[t+1] <- state
    }

  }
  EngageState <- data$Engage.State
  EngageState <- toupper(as.character(EngageState))
  data$Engage.State <- as.factor(EngageState)
  return(data)
}


read.Codes <- function(path="ConditionTimes.csv",startCol="Start.Time",endCol="End.Time",fmt="%H:%M:%S"){
if(grepl(".xls", path) || grepl(".xlsx", path)){
    data <- read_excel(path)
  }
  else {
    data <- read.csv(path,header = T)
  }
  
  print(as.character(data[,startCol]))
  data$Start.Time <- as.difftime(as.character(data[,startCol]),format = fmt,units = "secs")
  data$End.Time <- as.difftime(as.character(data[,endCol]),format = fmt,units = "secs")

  return(data)
}



read.ibutton <- function(path,header=T) {
  if(header){
    data <- read.csv(path,skip=14)
    
  }
  else {
    data <- read.csv(path,skip=14)
    
  }
  timeformat ="%m/%d/%y %I:%M:%S %p"
  data$Timestamp <- strptime(data$Date.Time,format=timeformat)
  data$Temperature <- data$Value
  return(subset(data,select=c("Timestamp","Temperature")))
}

read.reda <- function(file) {
  FIELDS <- c("Z","Y","X","Battery","Temperature","EDA","Tonic","Phasic","Events")
  #First lets read the file header info
  sample_rate <- as.numeric(scan(file, skip = 4, nlines = 1,sep = ":",strip.white = T,what = list("character","numeric"))[[2]]) # only 1 line after the skipped one
  st <- scan(file, skip = 5, nlines = 1,sep = " ",strip.white = T,what = "character") # only 1 line after the skipped one
  start_date <- st[3]
  start_time <- st[4]
  tz <- as.numeric(strsplit(st[5],":")[[1]][2])
  timestr <- sprintf("%s %s",start_date,start_time)
  timeformat ="%Y-%m-%d %H:%M:%S"
  start <- strptime(timestr,format=timeformat)
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  data <- as.data.frame(read.csv(file,header = F,sep = ",",skip=8,col.names=FIELDS))
  data <- na.omit(data)
  timestamps <- seq(from = start, by=dt , along.with = data$EDA)
  data$Timestamp <- timestamps
  #attr(data,"class") <- "eda"
  return(data)
}
read.eda <- function(file) {
  #First lets read the file header info
  sample_rate <- as.numeric(scan(file, skip = 4, nlines = 1,sep = ":",strip.white = T,what = list("character","numeric"))[[2]]) # only 1 line after the skipped one
  st <- scan(file, skip = 5, nlines = 1,sep = " ",strip.white = T,what = "character",) # only 1 line after the skipped one
  start_date <- st[3]
  start_time <- st[4]
  tz <- as.numeric(strsplit(st[5],":")[[1]][2])
  timestr <- sprintf("%s %s",start_date,start_time)
  timeformat ="%Y-%m-%d %H:%M:%S"
  start <- strptime(timestr,format=timeformat)
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  FIELDS <- scan(file, skip = 6, nlines = 1,sep = "|",strip.white = T,what = "character") # only 1 line after the skipped one
  FIELDS[1:6] <- c("Z","Y","X","Battery","Temperature","EDA")
  FIELDS[length(FIELDS)] <- strsplit(FIELDS[length(FIELDS)],split = " ")[[1]]
  data <- as.data.frame(read.csv(file,header = F,sep = ",",skip=8,col.names=FIELDS))
  data <- na.exclude(data)
  timestamps <- seq(from = start, by=dt , along.with = data$EDA)
  data$Timestamp <- timestamps
  data$Motion <- sqrt(data$X*data$X+data$Y*data$Y + data$Z*data$Z)
  #attr(data,"class") <- "eda"
  return(data)
}

read.actiwave <- function(file) {
  
  timeformat ="%m/%d/%Y %H:%M:%S"
  data <- as.data.frame(read.csv(file,header = T,sep = ","))
  data$Estimated_HR[data$Estimated_HR == 0] <- NA
  data <- na.omit(data)
  ts <- as.POSIXlt( strptime(data$Timestamp,format = timeformat)  )
  data$Timestamp <- ts
  data$IBI <- 1.0/(data$Estimated_HR/60.0)

  #attr(data,"class") <- "actiwave"
  
  return(data)
}

getFS <- function(d){
  fs <- 1.0/diff(d$Timestamp)[[1]]
  return(fs)
}


#/Users/OliverWS/Downloads/1418312035_192103
read.empatica <- function(path) {
  FIELDS <- c("Z","Y","X","Battery","Temperature","EDA")
  FILES <- c("ACC.csv","BVP.csv","EDA.csv","TEMP.csv","IBI.csv")
  #First lets read the file header info
  acc <- read.empatica.acc(file.path(path,"ACC.csv"))
  bvp <- read.empatica.bvp(file.path(path,"BVP.csv"))
  ibi <- read.empatica.ibi(file.path(path,"IBI.csv"))
  temp <- read.empatica.ibi(file.path(path,"IBI.csv"))
  eda <- read.empatica.ibi(file.path(path,"EDA.csv"))
  
  data <- list(ACC=acc,BVP=bvp,EDA=eda,TEMP=temp,IBI=ibi)
  attr(data,"class") <- "eda"
  return(data)
}

read.empatica.eda <- function(file){
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  sample_rate <- as.numeric(raw$V1[2])
  data <- data.frame(EDA=raw$V1[3:length(raw$V1)])
  start <- as.POSIXct(start_s, origin = "1970-01-01")
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  timestamps <- seq(from = start, by=dt , along.with = data$EDA) 
  data$Timestamp <- timestamps
  data
}

read.empatica.acc <- function(file){
  
}

read.empatica.bvp <- function(file){
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  sample_rate <- as.numeric(raw$V1[2])
  data <- data.frame(BVP=raw$V1[3:length(raw$V1)])
  start <- as.POSIXct(x = start_s, origin = "1970-01-01")
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  timestamps <- seq(from = start, by=dt , along.with = data$BVP) 
  data$Timestamp <- timestamps
  data
}

read.empatica.ibi <- function(file){
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  start <- as.POSIXct(x = start_s,origin = "1970-01-01")
  dt <- as.difftime(raw$V1[2:length(raw$V1)],units = "secs")
  timestamps <- start+dt
  ibi <- raw$V2[2:length(raw$V2)]
  data <- data.frame(Timestamp=timestamps, IBI=ibi)
  data
}

read.empatica.temp <- function(file) {
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  sample_rate <- as.numeric(raw$V1[2])
  temperatureF <- (9.0/5.0)*(raw$V1[3:length(raw$V1)] + 32.0)
  data <- data.frame(TEMP=temperatureF)
  start <- as.POSIXct(x = start_s, origin = "1970-01-01")
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  timestamps <- seq(from = start, by=dt , along.with = data$TEMP) 
  data$Timestamp <- timestamps
  data
}

as.simulateddyad <- function(p1,p2,cols=c("EDA"),norm=F) {
  p1.name <- deparse(substitute(p1))
  p2.name <- deparse(substitute(p2))
  
  p1.start <- 0
  p2.start <- 0
  p1.end <- length(p1$Timestamp)
  p2.end <- length(p2$Timestamp)
  
  dyad.start <- 0
  dyad.end <- min(c(p1.end,p2.end))
  if(dyad.end <= dyad.start){
    print("Error: p1 and p2 do not overlap")
    return()
  }
  print(paste("Dyad Overlap: ",dyad.start,"to",dyad.end))
  dyad.p1 <- p1[1:dyad.end,]
  dyad.p2 <- p2[1:dyad.end,]
  dyad <- data.frame(Timestamp=dyad.p1$Timestamp)
  for(c in cols) {
    cname.1 <- paste(p1.name,c,sep = ".")
    cname.2 <- paste(p2.name,c,sep = ".")
    if(norm){
      dyad[[cname.1]] <- normalize(dyad.p1[[c]])
      dyad[[cname.2]] <- normalize(dyad.p2[[c]])
    }
    else {
      dyad[[cname.1]] <- dyad.p1[[c]]
      dyad[[cname.2]] <- dyad.p2[[c]]
    }
    
  }
  
  return(dyad)
  
  
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

as.dyad <- function(p1,p2,cols=c("EDA"),norm=F,verbose=F) {
  p1.name <- deparse(substitute(p1))
  p2.name <- deparse(substitute(p2))
  
  p1.start <- min(p1$Timestamp)
  p2.start <- min(p2$Timestamp)
  p1.end <- max(p1$Timestamp)
  p2.end <- max(p2$Timestamp)
  
  dyad.start <- max(c(p1.start,p2.start))
  dyad.end <- min(c(p1.end,p2.end))
  if(dyad.end <= dyad.start){
    print("Error: p1 and p2 do not overlap")
    return()
  }
  if(verbose){
    print(paste("Dyad Overlap: ",dyad.start,"to",dyad.end))
    
  }
  dyad.p1 <- subset(p1,((Timestamp >= dyad.start) & (Timestamp < dyad.end)) )
  dyad.p2 <- subset(p2,((Timestamp >= dyad.start) & (Timestamp < dyad.end)) )
  l1 <- length(dyad.p1[[cols[1]]])
  l2 <- length(dyad.p2[[cols[1]]])
  if(l1 != l2){
    dyad.p1 <- dyad.p1[1:min(c(l1,l2)),]
    dyad.p2 <- dyad.p2[1:min(c(l1,l2)),]
  }
  dyad <- data.frame(Timestamp=dyad.p1$Timestamp)
  for(c in cols) {
    cname.1 <- paste(p1.name,c,sep = ".")
    cname.2 <- paste(p2.name,c,sep = ".")
    if(norm){
      dyad[[cname.1]] <- o.scale(dyad.p1[[c]])
      dyad[[cname.2]] <- o.scale(dyad.p2[[c]])
    }
    else {
      dyad[[cname.1]] <- dyad.p1[[c]]
      dyad[[cname.2]] <- dyad.p2[[c]]
    }

  }
  
  return(dyad)
  
  
}