library(readxl)
library(edfReader)
library(lubridate)

#' Read IBI Data from CSV
#'
#' This function reads interbeat interval (IBI) data from a CSV file.
#'
#' @param file File path to the CSV file.
#' @param startTime Starting time as a string (default "1970-01-01 0:0:0").
#' @param start Starting time as a POSIXlt object (default parsed from startTime).
#' @return A data frame with Timestamp and IBI columns.
#' @export
read.IBI <- function(file, startTime="1970-01-01 0:0:0",start=strptime(startTime,"%Y-%m-%d %H:%M:%S",tz = "")){
  raw <- read.csv(file,header = F)
  ibi.secs <- (raw$V1/1000.0)
  timestamps <- start + as.difftime(ibi.secs,units = "secs")
  IBI <- c(NA,diff(raw$V1))
  data <- data.frame(Timestamp=timestamps, IBI=IBI)
  data
}

#' Read IBI Time Series Data
#'
#' Reads IBI data and constructs a time series with computed heart rate (HR).
#'
#' @param file File path to the CSV file.
#' @param startTime Starting time string (default "1970-01-01 0:0:0").
#' @param start Starting time as POSIXlt (default parsed from startTime).
#' @return A data frame with Timestamp, IBI, and HR columns.
#' @export
read.IBI.TS <- function(file, startTime="1970-01-01 0:0:0",start=strptime(startTime,"%Y-%m-%d %H:%M:%S",tz = "")){
  data <- read.IBI(file,startTime = startTime,start=start)
  end <- max(data$Timestamp)
  dt <- as.difftime(as.character(1.0),format = "%OS")
  tstamps.secs <- seq(from = start, by=dt ,to = end)
  output <- data.frame(Timestamp=tstamps.secs, IBI=NA)
  for (n in seq_along(tstamps.secs)) {
    t = tstamps.secs[n]
    try(output$IBI[n] <- mean(data$IBI[which((data$Timestamp >= t) & (data$Timestamp < (t + dt)))],na.rm = T))
  }
  output$HR <- 60.0/(output$IBI/1000.0)
  return(output)
  
}


#' Read Engagement Codes
#'
#' Reads engagement codes from an Excel or CSV file.
#'
#' @param path File path to the codes file.
#' @return A data frame containing engagement codes.
#' @export
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


#' Read Condition Codes
#'
#' Reads condition codes from a CSV or Excel file.
#'
#' @param path File path (default "ConditionTimes.csv").
#' @param startCol Column name for start time.
#' @param endCol Column name for end time.
#' @param conditionCol Column name for condition.
#' @param fmt Format string for time (default "%H:%M:%S").
#' @param ... Additional parameters.
#' @return A data frame with condition codes.
#' @export
read.Codes <- function(path="ConditionTimes.csv",startCol="Start.Time",endCol="End.Time",conditionCol="Condition",fmt="%H:%M:%S",...){
if(grepl(".xls", path) || grepl(".xlsx", path)){
    data <- read_excel(path)
  }
  else {
    data <- read.csv(path,header = T,blank.lines.skip = T)
  }
  data <- na.omit(data)
  data$Start.Time <- as.difftime(as.character(data[,startCol]),format = fmt,units = "secs")
  data$End.Time <- as.difftime(as.character(data[,endCol]),format = fmt,units = "secs")
  data$Condition <-factor(data[,conditionCol], levels =data[,conditionCol][order(data[,startCol], data$Condition)], ordered=TRUE)
  data$Condition <- factor(data$Condition[,drop=TRUE])
  
  return(data)
}


#' Read RTCodes
#'
#' Reads RTCodes from a CSV or Excel file.
#'
#' @param path File path (default "ConditionTimes.csv").
#' @param startCol Column name for start time.
#' @param endCol Column name for end time.
#' @param fmt Time format (default "%y-%m-%d %H:%M:%S").
#' @return A data frame with RTCodes.
#' @export
read.RTCodes <- function(path="ConditionTimes.csv",startCol="Start.Time",endCol="End.Time",fmt="%y-%m-%d %H:%M:%S"){
  if(grepl(".xls", path) || grepl(".xlsx", path)){
    data <- read_excel(path)
  }
  else {
    data <- read.csv(path,header = T,blank.lines.skip = T)
  }
  
  data <- na.omit(data)
  data$Start.Time <- as.character.POSIXt(as.character(data[,startCol]),format = fmt)
  data$End.Time <- as.character.POSIXt(as.character(data[,endCol]),format = fmt)
  data$Condition <-factor(data$Condition, levels =data$Condition[order(data[,startCol], data$Condition)], ordered=TRUE)
  data$Condition <- factor(data$Condition[,drop=TRUE])
  
  return(data)
}




#' Read iButton Data
#'
#' Reads data from an iButton CSV file skipping header lines.
#'
#' @param path File path to the CSV file.
#' @param header Logical indicating if file has header (default TRUE).
#' @return A data frame with Timestamp and Temperature columns.
#' @export
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

#' Read Reda Data
#'
#' Reads and processes REDA data from a CSV file.
#'
#' @param file File path to the REDA data file.
#' @return A data frame containing the processed REDA data.
#' @export
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


#' Read EDA Data
#'
#' Reads EDA (electrodermal activity) data from a CSV file.
#'
#' @param file File path to the EDA data file.
#' @return A data frame with EDA data, Timestamp, and Motion computed.
#' @export
read.eda <- function(file) {
  #First lets read the file header info
  sample_rate <- as.numeric(scan(file, skip = 4, nlines = 1,sep = ":",strip.white = T, quiet=T, what = list("character","numeric"))[[2]]) # only 1 line after the skipped one
  st <- scan(file, skip = 5, nlines = 1,sep = " ",strip.white = T, quiet=T, what = "character",) # only 1 line after the skipped one
  start_date <- st[3]
  start_time <- st[4]
  tz <- as.numeric(strsplit(st[5],":")[[1]][2])
  timestr <- sprintf("%s %s",start_date,start_time)
  timeformat ="%Y-%m-%d %H:%M:%S"
  start <- strptime(timestr,format=timeformat)
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  FIELDS <- scan(file, skip = 6, nlines = 1,sep = "|",strip.white = T, quiet=T, what = "character") # only 1 line after the skipped one
  FIELDS[1:6] <- c("Z","Y","X","Battery","Temperature","EDA")
  FIELDS[length(FIELDS)] <- strsplit(FIELDS[length(FIELDS)],split = " ")[[1]]
  data <- as.data.frame(read.csv(file,header = F,sep = ",",skip=8,col.names=FIELDS,))
  data <- na.exclude(data)
  timestamps <- seq(from = start, by=dt , along.with = data$EDA)
  data$Timestamp <- timestamps
  data$Motion <- sqrt(data$X*data$X+data$Y*data$Y + data$Z*data$Z)
  #attr(data,"class") <- "eda"
  return(data)
}

#' Save EDA Data to File
#'
#' Writes EDA data to a file with a header.
#'
#' @param data Data frame containing EDA data.
#' @param filename Output file name.
#' @export
save.eda <- function(data, filename) {
  # Log File Created by Q Analytics - (c) 2011 Affectiva Inc.
  # File Version: 1.01
  # Firmware Version: 1.81
  # UUID: AQL181320DU
  # Sampling Rate: 32
  # Start Time: 2015-05-12 16:31:47 Offset:-04
  # Z-axis | Y-axis | X-axis | Battery | °Celsius | EDA(uS) 
  # ---------------------------------------------------------
  
  sampleRate = getFS(data)
  startTime = min(data$Timestamp)
  timeString = strftime(startTime, format="%Y-%m-%d %H:%M:%S Offset:%z")
  timeString = strtrim(timeString,str_length(timeString)-2)
  header = c(
    "Log File Created by CBS Toolkit",
    "File Version: 1.01",
    "Firmware Version: 1.81",
    "UUID: AQL0000000",
    paste("Sampling Rate:", round(sampleRate),sep=" "),
    paste("Start Time:",timeString, sep=" "),
    "Z-axis | Y-axis | X-axis | Battery | °Celsius | EDA(uS)",
    "---------------------------------------------------------"
  )
  writeLines(header,con=filename,sep="\r\n")
  
  data.sub <- subset(data,select = c("Z","Y","X","Battery","Temperature","EDA"))
  write.table(data.sub,file = filename,append = T,sep = ",",row.names = F, col.names = F, quote = F,eol = "\r\n",na = "-1.0")
  
  
}
write.eda <- save.eda

#' Read CSV Data with Timestamp Parsing
#'
#' Reads a CSV file and parses the Timestamp column using ymd_hms.
#'
#' @param file File path to the CSV file.
#' @return A data frame with a parsed Timestamp column.
#' @export
read.csvdata <- function(file) {
  
  timeformat ="%Y-%m-%d %H:%M:%OS"
  data <- as.data.frame(read.csv(file,header = T,sep = ","))
  data$Timestamp <- ymd_hms(data$Timestamp,tz = "EST5EDT")
  return(data)
}


#' Read Biosync Data
#'
#' Reads biosync data from a CSV file and fixes the Timestamp format.
#'
#' @param file File path to the CSV file.
#' @return A data frame with a corrected Timestamp column.
#' @export
read.biosync <- function(file) {
  
  timeformat ="%Y-%m-%d %H:%M:%OS"
  data <- as.data.frame(read.csv(file,header = T,sep = ","))
  data$Timestamp <- gsub("(.*)\\:", "\\1.", data$Timestamp)
  data$Timestamp <- ymd_hms(data$Timestamp,tz = "EST5EDT")
  return(data)
}

#' Read EDF Data
#'
#' Reads EDF (European Data Format) data using edfReader.
#'
#' @param file File path to the EDF file.
#' @return A list of data frames, one per signal.
#' @export
read.edf <- function(file){
  header = edfReader::readEdfHeader(file)
  data = edfReader::readEdfSignals(header,signals = "Ordinary")
  
  output <- list()
  for(signal in data){
    safe_name <- make.names(signal$label)
    start <- strptime(signal$startTime, format="%Y-%m-%d %H:%M:%S")
    sample_rate <- signal$sRate
    dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
    timestamps <- seq(from = start, by=dt , along.with = signal$signal) 
    df <- data.frame(Timestamp=timestamps)
    df[[safe_name]] <- signal$signal
    output[[safe_name]] <- df
  }
  return(output)
}
#' Read Actiwave Data
#'
#' Reads Actiwave data from a CSV file and computes IBI.
#'
#' @param file File path to the Actiwave CSV file.
#' @return A data frame with Timestamp, HR, and IBI columns.
#' @export
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

#' Get Sampling Frequency
#'
#' Computes the sampling frequency from the Timestamp column of a data frame.
#'
#' @param d Data frame with a Timestamp column.
#' @return The sampling frequency.
#' @export
getFS <- function(d){
  fs <- 1.0/diff(d$Timestamp)[[1]]
  if(fs < 1.0){
    return(fs)
  }
  else {
    return(round(fs))
  }
}


#' Downsample EDA Data
#'
#' Downsamples EDA data to a target sampling rate.
#'
#' @param d Data frame with EDA data.
#' @param targetRate Desired sampling rate (default 1.0).
#' @return A downsampled data frame.
#' @export
downsample.eda <- function(d,targetRate=1.0){
  current.fs <- getFS(d)
  q <- round(current.fs/targetRate)
  output.TS <- d$Timestamp[seq(from=1,to=dim(d)[1],by=q)]
  output <- data.frame(Timestamp=output.TS)
  for(i in 1:dim(d)[2]){
    col <- colnames(d)[i]
    if(col == "Timestamp"){
      #Skip this column
    }
    else {
      x <- d[,i]
      output[,col] <- decimate(x,q)
    }
  }
  return(output)
}



#/Users/OliverWS/Downloads/1418312035_192103
#' Read Empatica Data
#'
#' Reads Empatica device data from a specified directory.
#'
#' @param path Directory path containing Empatica CSV files.
#' @return A list of data frames for different signals (ACC, BVP, EDA, TEMP, IBI).
#' @export
read.empatica <- function(path) {
  FIELDS <- c("Z","Y","X","Battery","Temperature","EDA")
  FILES <- c("ACC.csv","BVP.csv","EDA.csv","TEMP.csv","IBI.csv")
  #First lets read the file header info
  acc <- read.empatica.acc(file.path(path,"ACC.csv"))
  bvp <- read.empatica.bvp(file.path(path,"BVP.csv"))
  ibi <- read.empatica.ibi(file.path(path,"IBI.csv"))
  temp <- read.empatica.temp(file.path(path,"TEMP.csv"))
  eda <- read.empatica.eda(file.path(path,"EDA.csv"))
  
  data <- list(ACC=acc,BVP=bvp,EDA=eda,TEMP=temp,IBI=ibi)
  attr(data,"class") <- "eda"
  return(data)
}

#' Read E4 Data
#'
#' Reads E4 device data from a specified directory.
#'
#' @param path Directory path containing E4 CSV files.
#' @return A data frame with processed E4 data.
#' @export
read.e4 <- function(path) {
  FIELDS <- c("Z","Y","X","Battery","Temperature","EDA")
  FILES <- c("ACC.csv","BVP.csv","EDA.csv","TEMP.csv","IBI.csv")
  #First lets read the file header info
  acc <- read.empatica.acc(file.path(path,"ACC.csv"))
  bvp <- read.empatica.bvp(file.path(path,"BVP.csv"))
  try(ibi <- read.empatica.ibi(file.path(path,"IBI.csv")))
  temp <- read.empatica.temp(file.path(path,"TEMP.csv"))
  eda <- read.empatica.eda(file.path(path,"EDA.csv"))
  hr <- read.empatica.hr(file.path(path,"HR.csv"))
  df <- data.frame(Timestamp=bvp$Timestamp)
  df$EDA <- interp(eda$EDA,q = 8,n = 5)
  df$BVP <- bvp$BVP
  df$X <- interp(acc$X,q = 2,n = 5)
  df$Y <- interp(acc$Y,q = 2,n = 5)
  df$Z <- interp(acc$Z,q = 2,n = 5)
  df$Temp <- interp(temp$Temp,q = 8,n = 5)
  df$HR <- interp(temp$HR,q = 64,n = 5)
  return(df)
}


#' Read Empatica EDA Data
#'
#' Reads EDA data from an Empatica EDA CSV file.
#'
#' @param file File path to the EDA CSV file.
#' @return A data frame with Timestamp and EDA data.
#' @export
read.empatica.eda <- function(file){
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  sample_rate <- as.numeric(raw$V1[2])
  data <- data.frame(Timestamp=NA, EDA=raw$V1[3:length(raw$V1)])
  start <- as.POSIXct(start_s, origin = "1970-01-01")
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  timestamps <- seq(from = start, by=dt , along.with = data$EDA) 
  data$Timestamp <- timestamps
  data
}

#' Read Empatica Accelerometer Data
#'
#' Reads accelerometer data from an Empatica ACC CSV file.
#'
#' @param file File path to the ACC CSV file.
#' @return A data frame with X, Y, Z, and Timestamp columns.
#' @export
read.empatica.acc <- function(file){
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  sample_rate <- as.numeric(raw$V1[2])
  data <- data.frame(X=raw$V1[3:length(raw$V1)]/64.0,Y=raw$V2[3:length(raw$V2)]/64.0,Z=raw$V3[3:length(raw$V3)]/64.0)
  start <- as.POSIXct(x = start_s, origin = "1970-01-01")
  dt <- as.difftime(as.character(1.0/sample_rate),format = "%OS")
  timestamps <- seq(from = start, by=dt , along.with = data$X) 
  data$Timestamp <- timestamps
  data
}

#' Read Empatica BVP Data
#'
#' Reads blood volume pulse (BVP) data from an Empatica BVP CSV file.
#'
#' @param file File path to the BVP CSV file.
#' @return A data frame with BVP and Timestamp columns.
#' @export
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
read.empatica.hr <- function(file){
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  start <- as.POSIXct(x = start_s,origin = "1970-01-01")
  dt <- as.difftime(seq_along(raw$V1[2:length(raw$V1)]),units = "secs")
  timestamps <- start+dt
  HR <- as.numeric(raw$V1[2:length(raw$V1)])
  data <- data.frame(Timestamp=timestamps, HR=HR)
  data
}

read.empatica.ibi <- function(file){
  raw <- read.csv(file,header = F)
  start_s <- raw$V1[1]
  start <- as.POSIXct(x = start_s,origin = "1970-01-01")
  dt <- as.difftime(raw$V1[2:length(raw$V1)],units = "secs")
  timestamps <- start+dt
  ibi <- as.double(as.character(raw$V2[2:length(raw$V2)]) )
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

#' Create Simulated Dyad Data
#'
#' Simulates a dyad using two participant data frames.
#'
#' @param p1 Data frame for participant 1.
#' @param p2 Data frame for participant 2.
#' @param cols Columns to include (default "EDA").
#' @param norm Logical indicating whether to normalize data (default FALSE).
#' @param verbose Logical for verbose output (default FALSE).
#' @return A data frame representing the simulated dyad.
#' @export
as.simulateddyad <- function(p1,p2,cols=c("EDA"),norm=F,verbose=F) {
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
  if(verbose){
    print(paste("Dyad Overlap: ",dyad.start,"to",dyad.end))
  }
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

#' Scale Vector to [0,1]
#'
#' Scales a numeric vector to the range [0,1].
#'
#' @param x Numeric vector.
#' @param use.sd Logical, if TRUE uses standard deviation for scaling (default FALSE).
#' @return Scaled numeric vector.
#' @export
o.scale <- function(x,use.sd=F){
  if(use.sd){
    sx <- (x - min(x, na.rm=T))/sd(x,na.rm = T)
    
  }
  else {
    sx <- (x - min(x, na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))    
  }
  return(sx)
  
}




#' Read Dyad Data from File
#'
#' Reads dyad data from a file and constructs a data frame with Timestamps assigned based on the sample rate.
#'
#' @param f File path to the dyad data file.
#' @param sample_rate Sampling rate (default 32.0).
#' @param p1.name Optional name for participant 1 column (if NULL, first column name is used).
#' @param p2.name Optional name for participant 2 column (if NULL, second column name is used).
#' @param start Starting time as POSIXlt (default "2016-01-01 0:0:0").
#' @param readFunction Function to read the file (default read.csv).
#' @return A data frame with Timestamp and dyad data.
#' @export
read.dyad <- function(f,sample_rate=32.0,p1.name=NULL,p2.name=NULL,start=strptime("2016-01-01 0:0:0","%Y-%m-%d %H:%M:%S"), readFunction=read.csv){
  data <- as.data.frame(readFunction(f))
  if(is.null(p1.name)){
    p1.name <- colnames(data)[1]
  }
  if(is.null(p2.name)){
    p2.name <- colnames(data)[2]
  }
  
  sample_period=1.0/sample_rate
  nSamples = dim(data)[1]
  dt <- as.difftime(sample_period,units = "secs")
  timestamps <- seq(from = start, by=dt , length.out = nSamples)
  outputData <- data.frame(Timestamp = timestamps)
  outputData[,2] <- data[,1]
  outputData[,3] <- data[,2]
  colnames(outputData) <- c("Timestamp",p1.name,p2.name)
  return(na.omit(outputData))
}

#' Convert Two Participant Data into a Dyad
#'
#' Aligns and merges two participant data frames into a dyadic data frame based on overlapping timestamps.
#'
#' @param p1 Data frame for participant 1.
#' @param p2 Data frame for participant 2.
#' @param cols Columns to include (default "EDA").
#' @param norm Logical, if TRUE normalize the data (default FALSE).
#' @param verbose Logical for verbose output (default FALSE).
#' @param na.interpolate Logical, if TRUE interpolate NA values (default TRUE).
#' @param interpolation.method Function to interpolate NAs (default na.spline).
#' @return A dyadic data frame with aligned timestamps.
#' @export
as.dyad <- function(p1,p2,cols=c("EDA"),norm=F,verbose=F,na.interpolate=T,interpolation.method=na.spline) {
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
    
    if(na.interpolate){
      dyad[[cname.1]] <- interpolation.method(dyad[[cname.1]])
      dyad[[cname.2]] <- interpolation.method(dyad[[cname.2]])
      
    }

  }
  
  return(dyad)
  
  
}