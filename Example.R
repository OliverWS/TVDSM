source("~/git/CBSL.R/R/utils.R")
source("~/git/CBSL.R/R/read.R")
source("~/git/CBSL.R/R/plot.R")
source("~/git/CBSL.R/R/physio.R")

setwd("~/git/CBSL.R/")
child <- read.eda(file = "Child_cleaned.eda" )
therapist <- read.eda(file="Therapist_cleaned.eda")

myDyad <- as.dyad(child,therapist)


#mdl <- statespace(myDyad$child.EDA,myDyad$therapist.EDA,type = 2)

computeStateSpace <- function(dyad) {
  mdl <- statespace(decimate(dyad$child.EDA,32),decimate(dyad$therapist.EDA, 32),type = 3)
  return(mdl)
}

w <- o.window.list(myDyad,window_size = 32*60*2,FUN = computeStateSpace,na.rm = T)




