### Imports data from CSV files in our couponbonds list format

## NOTE: This function is experimental and we only provide it
##       because we have been asked several times how we do the
##       data import.
##       A detailed description of our the format is given in
##       the paper, so a custom import function depending on the
##       the data source should be easy to code.

## RESTRICTIONS: only annual coupons, no bond should mature within
##               sample range

dyncouponbonds <- function(datafiles,   # vector with CSV files c("COUNTRY S.csv", "COUNTRY AC.csv", "COUNTRY CP.csv")
                           gname,       # country name
                           export=FALSE # write imported data to RData file
                           ) {

  rawdata <- read.csv(datafiles[1], dec=".",sep=",",colClasses = "character")
  rawdata_AC <- read.csv(datafiles[2], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
  rawdata_CP <- read.csv(datafiles[3], dec=".",sep=",",colClasses = "character", header=FALSE, check.names=FALSE)
 
  DATES <- as.character(as.Date(rawdata_AC[,1],format="%d%m%Y"))
  DATES <- paste(substring(DATES,9,10),substring(DATES,6,7),substring(DATES,1,4),sep="")


  ## dates where accrued interest is acutally paid (only works if all list elements are of lenght one)
  AC_dates <- list()
  for(i in 1:dim(rawdata)[1]){
    AC_dates[[i]] <- rawdata_AC[which(diff(as.numeric(rawdata_AC[,-1][,i]))<0)+1,1]
    if (length(AC_dates[[i]]) == 0){    # empty string for bonds with no coupon payment in dataset
      AC_dates[[i]] <- ""
    }
    
  }  
  names(AC_dates) <- rawdata[,1]
  dslist <- list()

  
  for(i in 1:length(DATES)) {

    TODAY <- DATES[i]
    AC <- rawdata_AC[,-1][i,]
    CP <- rawdata_CP[,-1][i,]
    
    datastreamlist <- function(rawdata, TODAY, AC, CP, AC_dates) {
      
      data <- list()
      data$ISIN <- rawdata[,"ISIN"]
      data$MATURITYDATE <- as.Date(rawdata[,"MATURITYDATE"],format="%d%m%Y")
      data$ISSUEDATE <- as.Date(rawdata[,"ISSUEDATE"],format="%d%m%Y")
      data$COUPONRATE <- as.numeric(rawdata[,"COUPONRATE"])
      data$PRICE <- as.numeric(CP)
      data$ACCRUED <- as.numeric(AC)
      data$CASHFLOWS <- list()
      data$TODAY <- as.Date(TODAY,format="%d%m%Y")
      
      NEXTCOUPON <- ifelse(as.Date(paste(rawdata[,"COUPONDATE"],substring(TODAY,5,8),sep=""),format="%d%m%Y") > data$TODAY,
                           paste(rawdata[,"COUPONDATE"],substring(TODAY,5,8),sep=""),
                           paste(rawdata[,"COUPONDATE"],as.character(as.numeric(substring(TODAY,5,8))+1),sep=""))
      NCOUPON <- as.numeric(substring(rawdata[,"MATURITYDATE"],5,8)) - as.numeric(substring(NEXTCOUPON,5,8))

      ## cash flows ISIN
      data$CASHFLOWS$ISIN <- vector()
      for(i in 1:length(NCOUPON)){
        data$CASHFLOWS$ISIN <-  c(data$CASHFLOWS$ISIN,rep(data$ISIN[i],NCOUPON[i]+1))
      }
      
      ## cash flows
      data$CASHFLOWS$CF <- vector()
      for(i in 1:length(NCOUPON)){
        data$CASHFLOWS$CF <- c(data$CASHFLOWS$CF,c(rep(data$COUPONRATE[i]*100,NCOUPON[i]),100+data$COUPONRATE[i]*100))
      }
      
      ## cash flow dates
      data$CASHFLOWS$DATE <- vector()
      for(i in 1:length(NCOUPON)){
        data$CASHFLOWS$DATE <- c(data$CASHFLOWS$DATE,paste(rawdata[i,"COUPONDATE"],as.numeric(substring(NEXTCOUPON[i],5,8)) + seq(0,NCOUPON[i]),sep=""))
      }
      
      ## correct for AC dates
      ## get indices where new ISINs start
      cp_ind <- c(1,1+cumsum(1+NCOUPON))
      cp_ind <- cp_ind[-length(cp_ind)]
      
      for(j in 1:length(cp_ind)) {
        if(length(AC_dates[[j]]) > 1){
          if(as.numeric(substring(data$CASHFLOWS$DATE[cp_ind[j]],5,8)) <= as.numeric(substring(AC_dates[[j]],5,8))) data$CASHFLOWS$DATE[cp_ind[j]] = AC_dates[[j]]
        }
      }
      
      ## throw out cashflows at coupon date
      to_ind <- which(as.Date(data$CASHFLOWS$DATE,format="%d%m%Y") <= as.Date(TODAY,format="%d%m%Y"))
      if(length(to_ind)>0) {
        data$CASHFLOWS$ISIN <- data$CASHFLOWS$ISIN[-to_ind]
        data$CASHFLOWS$DATE <- data$CASHFLOWS$DATE[-to_ind]
        data$CASHFLOWS$CF <- data$CASHFLOWS$CF[-to_ind]
      }
      data$CASHFLOWS$DATE <- as.Date(data$CASHFLOWS$DATE,format="%d%m%Y")
      data
    }
    
    dslist[[i]] <- list(datastreamlist(rawdata, TODAY, AC, CP, AC_dates))
    names(dslist[[i]]) <- gname
    class(dslist[[i]]) <- "couponbonds"
  }
  class(dslist) <- "dyncouponbonds"
  save(dslist,file=paste(gname,".RData",sep=""))
  
  dslist
}

print.dyncouponbonds <- function(x, ...) {
  cat("This is a dynamic dataset of coupon bonds.\n")
  cat(paste("There are",length(x), "observations between",x[[1]][[1]]$TODAY, "and",x[[length(x)]][[1]]$TODAY,".\n"))
  }

print.couponbonds <- function(x, ...) {
  cat("This is a dataset of coupon bonds for:\n")
  cat(names(x),",","\n")
  cat(paste("observed at ", x[[1]]$TODAY,".","\n",sep=""))
  ## TODO: more info
  }
