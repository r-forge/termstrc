###################################################################
#               Import Datastream data                            #
###################################################################

# RESTRICTIONS:
# - no bond is allowed to mature within the sample range
# - only annual coupons
# - data set should not be longer than a year, otherwise problems
#   with AC dates

import_ds <- function(datafiles,gname,export=TRUE) {

#################
#rm(list=ls(all=TRUE))
#datafiles <- c("DataGermany S.csv", "DataGermany AC.csv", "DataGermany CP.csv", "DataGermany VO.csv", "DataGermany SETT.csv", "DataGermany RY.csv")
#gname <- "GERMANY"
#################


rawdata <- read.csv(datafiles[1], dec=".",sep=",",colClasses = "character")
rawdata_AC <- read.csv(datafiles[2], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
rawdata_CP <- read.csv(datafiles[3], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
rawdata_VO <- read.csv(datafiles[4], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
rawdata_SETT <- read.csv(datafiles[5], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
rawdata_RY <- read.csv(datafiles[6], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
rawdata_IBOX <- read.csv(datafiles[7], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
rawdata_IBXA <- read.csv(datafiles[8], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
rawdata_IBXB <- read.csv(datafiles[9], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
 
DATES <- as.character(as.Date(rawdata_AC[,1],format="%d%m%Y"))
DATES <- paste(substring(DATES,9,10),substring(DATES,6,7),substring(DATES,1,4),sep="")


# dates where accrued interst is acutally paid (only works if all list elements are of lenght one)
AC_dates <- list()
for(i in 1:dim(rawdata)[1]){
  AC_dates[[i]] <- rawdata_AC[which(diff(as.numeric(rawdata_AC[,-1][,i]))<0)+1,1]
print(AC_dates[[i]])
}  
names(AC_dates) <- rawdata[,1]

dslist <- list()

for(i in 1:length(DATES)) {
 # for(i in 44:45){
  TODAY <- DATES[i]
  AC <- rawdata_AC[,-1][i,]
  CP <- rawdata_CP[,-1][i,]
  VO <- rawdata_VO[,-1][i,]
  RY <- rawdata_RY[,-1][i,]
  IBOX <- rawdata_IBOX[,-1][i,]
  IBXA <- rawdata_IBXA[,-1][i,]
  IBXB <- rawdata_IBXB[,-1][i,]
  
datastreamlist <- function(rawdata, TODAY, AC, CP, VO, RY, AC_dates,IBOX,IBXA,IBXB) {
  data <- list()
  data$ISIN <- rawdata[,"ISIN"]
  data$MATURITYDATE <- as.Date(rawdata[,"MATURITYDATE"],format="%d%m%Y")
  data$ISSUEDATE <- as.Date(rawdata[,"ISSUEDATE"],format="%d%m%Y")
  data$COUPONRATE <- as.numeric(rawdata[,"COUPONRATE"])
  data$PRICE <- as.numeric(CP)
  data$VOLUME <- as.numeric(VO)
  data$RY <- as.numeric(RY)
  data$IBOX <- as.numeric(IBOX)
  data$IBXA <- as.numeric(IBXA)
  data$IBXB <- as.numeric(IBXB)
  data$ACCRUED <- as.numeric(AC)
  data$CASHFLOWS <- list()
  data$TODAY <- as.Date(TODAY,format="%d%m%Y")

  NEXTCOUPON <- ifelse(as.Date(paste(rawdata[,"COUPONDATE"],substring(TODAY,5,8),sep=""),format="%d%m%Y") > data$TODAY,
                       paste(rawdata[,"COUPONDATE"],substring(TODAY,5,8),sep=""),
                       paste(rawdata[,"COUPONDATE"],as.character(as.numeric(substring(TODAY,5,8))+1),sep=""))
  NCOUPON <- as.numeric(substring(rawdata[,"MATURITYDATE"],5,8)) - as.numeric(substring(NEXTCOUPON,5,8))

  # cash flows ISIN
  data$CASHFLOWS$ISIN <- vector()
  for(i in 1:length(NCOUPON)){
    data$CASHFLOWS$ISIN <-  c(data$CASHFLOWS$ISIN,rep(data$ISIN[i],NCOUPON[i]+1))

    }

  # cash flows
  data$CASHFLOWS$CF <- vector()
  for(i in 1:length(NCOUPON)){
    data$CASHFLOWS$CF <- c(data$CASHFLOWS$CF,c(rep(data$COUPONRATE[i]*100,NCOUPON[i]),100+data$COUPONRATE[i]*100))
    }

  # cash flow dates
  data$CASHFLOWS$DATE <- vector()
  for(i in 1:length(NCOUPON)){
    data$CASHFLOWS$DATE <- c(data$CASHFLOWS$DATE,paste(rawdata[i,"COUPONDATE"],as.numeric(substring(NEXTCOUPON[i],5,8)) + seq(0,NCOUPON[i]),sep=""))
    }

  # correct for AC dates
  # get indices where new ISINs start
  cp_ind <- c(1,1+cumsum(1+NCOUPON))
  cp_ind <- cp_ind[-length(cp_ind)]


  # WARNING: only works if there is no case at turn of the year
  for(j in 1:length(cp_ind)) {
    if(as.numeric(substring(data$CASHFLOWS$DATE[cp_ind[j]],5,8)) <= as.numeric(substring(AC_dates[[j]],5,8))) data$CASHFLOWS$DATE[cp_ind[j]] = AC_dates[[j]]
  }
  
  # throw out cashflows at coupon date
  to_ind <- which(as.Date(data$CASHFLOWS$DATE,format="%d%m%Y") <= as.Date(TODAY,format="%d%m%Y"))
  if(length(to_ind)>0) {
    data$CASHFLOWS$ISIN <- data$CASHFLOWS$ISIN[-to_ind]
    data$CASHFLOWS$DATE <- data$CASHFLOWS$DATE[-to_ind]
    data$CASHFLOWS$CF <- data$CASHFLOWS$CF[-to_ind]
  }

  
  data$CASHFLOWS$DATE <- as.Date(data$CASHFLOWS$DATE,format="%d%m%Y")

  data
}

  
  dslist[[i]] <- datastreamlist(rawdata, TODAY, AC, CP, VO, RY, AC_dates)  
}


  
names(dslist) <- rep(gname,length(dslist))

save(dslist,file=paste(gname,".RData",sep=""))

  dslist

}
