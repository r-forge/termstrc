import_ds <- function(datafilesn,gname,export=FALSE){
#datafiles = c("france S.csv", "france AC.csv", "france CP.csv")
 
rawdata <- read.csv(datafiles[1], dec=".",sep=",",colClasses = "character")
rawdata_AC <- read.csv(datafiles[2], dec=".",sep=",",colClasses = "character", header=TRUE, check.names=FALSE)
rawdata_CP <- read.csv(datafiles[3], dec=".",sep=",",colClasses = "character", header=TRUE, check.names=FALSE)

DATES <- as.character(as.Date(names(rawdata_AC)[-(1:2)], format="%d.%m.%y"))
DATES <- paste(substring(DATES,9,10),substring(DATES,6,7),substring(DATES,1,4),sep="")

AC <- rawdata[,"ACCRUED"]
CP <- rawdata[,"PRICE"]

dslist <- list()

for(i in 1:length(DATES)){
  TODAY <- DATES[i]
  AC <- rawdata_AC[-(1:2)][,i]
  CP <- rawdata_CP[-(1:2)][,i]
  
datastreamlist <- function(rawdata, TODAY, AC, CP){
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
 
  data$CASHFLOWS$DATE <- as.Date(data$CASHFLOWS$DATE,format="%d%m%Y")
  
data
}
  dslist[[i]] <- datastreamlist(rawdata, TODAY, AC, CP)
}
names(dslist) <- rep(gname,length(dslist))
if(is.character(export)) save(dslist,file=paste(export,".RData",sep=""))
dslist

}  

datafiles = c("france S.csv", "france AC.csv", "france CP.csv")
test <- import_ds(datafiles,gname="FRANCE",export="France_test")


 dyn_rm_bond <- function(bdata, ISIN) 
{
    for (i in seq(length(bdata))) {
     cf_isin_index <- which(bdata[[i]]$CASHFLOWS$ISIN %in% ISIN)
     isin_index <- which(bdata[[i]]$ISIN %in% ISIN)
     bdata[[i]]$ISIN <- bdata[[i]]$ISIN[-isin_index]
     bdata[[i]]$MATURITYDATE <- bdata[[i]]$MATURITYDATE[-isin_index]
     bdata[[i]]$STARTDATE <- bdata[[i]]$STARTDATE[-isin_index]
     bdata[[i]]$COUPONRATE <- bdata[[i]]$COUPONRATE[-isin_index]
     bdata[[i]]$PRICE <- bdata[[i]]$PRICE[-isin_index]
     bdata[[i]]$ACCRUED <- bdata[[i]]$ACCRUED[-isin_index]
     bdata[[i]]$CASHFLOWS$ISIN <- bdata[[i]]$CASHFLOWS$ISIN[-cf_isin_index]
     bdata[[i]]$CASHFLOWS$CF <- bdata[[i]]$CASHFLOWS$CF[-cf_isin_index]
     bdata[[i]]$CASHFLOWS$DATE <- bdata[[i]]$CASHFLOWS$DATE[-cf_isin_index]
   }
    bdata
}




dynbonddata <- dyn_rm_bond(test,c("FR0000571044","FR0000571085","FR0000570780","FR0107369672","FR0106589437","FR0108354806","FR0109970386"))


matrange <- c(0,20)
method <- "Svensson"
fit <- "prices"
weights <- "duration"

startparam <-  b <- matrix(c(0.04835332, -0.008060085, 0.3879382, 2.322683, -0.4039001, 2.353533),nrow=1,ncol=6,byrow=TRUE)


dyntermstrc <- function(dynbonddata,matrange,method,fit,weights,startparam,...) {
x <- list()
for (i in seq(length(dynbonddata))) {


if(i>1) b <- matrix(x[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=6,byrow=TRUE) else b <- startparam

group <- names(dynbonddata)[i]

rownames(b) <- group
colnames(b) <- c("beta0","beta1","beta2","tau1","beta3","tau2")

bonddata <- list()
bonddata[[group]] <- dynbonddata[[i]]
x[[i]] <- nelson_estim(group, bonddata= bonddata, matrange, 
                  method, fit, weights, startparam=b,...)

class(x[[i]]) <- "nelson"
}

class(x) <- "dyntermstrc"
x
}

myres  <- dyntermstrc(dynbonddata,matrange,method,fit,weights,startparam)



opt_result <- x[[1]]$opt_result$FRANCE$par
for (i in 2:length(govbondsts)){
  opt_result <- rbind(opt_result, x[[i]]$opt_result$FRANCE$par)
}

