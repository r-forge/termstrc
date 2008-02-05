
importdatastream <- function(datafiles){
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
dslist
}  

library("termstrc")

datafiles = c("france S.csv", "france AC.csv", "france CP.csv")
govbondsts <- importdatastream(datafiles)

x <- list()
for (i in 1:length(govbondsts)){
govbonds <- list()
govbonds$FRANCE <- govbondsts[[i]]

group <- c("FRANCE")
bonddata <- govbonds
matrange <- c(1,20)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "duration"
control <- list(eval.max=100000, iter.max=500)

if(i>1) b <- matrix(x[[i-1]]$opt_result$FRANCE$par,nrow=1,ncol=4,byrow=TRUE) else b <- matrix(c(0.1,0.1,0.1, 1),nrow=1,ncol=4,byrow=TRUE)

			
rownames(b) <- group
colnames(b) <- c("beta0","beta1","beta2","tau1")

x[[i]] <- nelson_estim(group, bonddata, matrange, 
                  method, fit, weights, startparam=b,control)
}

opt_result <- x[[1]]$opt_result$FRANCE$par
for (i in 2:length(govbondsts)){
  opt_result <- rbind(opt_result, x[[i]]$opt_result$FRANCE$par)
}

par(mfrow = c(2,2))
plot(opt_result[,1],type="l",ylab="beta_0", col=1, lwd=2)
grid()
plot(opt_result[,2],type="l",ylab="beta_1", col=2, lwd=2)
grid()
plot(opt_result[,3],type="l",ylab="beta_2", col=3, lwd=2)
grid()
plot(opt_result[,4],type="l",ylab="tau_1", col=4, lwd=2)
grid()


plot(seq(0,20,0.01), spotrates(method="Nelson/Siegel",opt_result[1,],seq(0,20,0.01)),type="l",ylim=c(0.03,0.07))

for (i in 2:nrow(opt_result)){
lines(seq(0,20,0.01), spotrates(method="Nelson/Siegel",opt_result[i,],seq(0,20,0.01)))
  }

library("rgl")

X <- seq(1.5,30,0.1)
Y <- seq(nrow(opt_result))
Z <- matrix(spotrates(method="Nelson/Siegel",opt_result[1,],X),nrow=1)
for (i in 2:nrow(opt_result)) {
  Z <- rbind(Z, spotrates(method="Nelson/Siegel",opt_result[i,],X))}

persp(X,Y,t(Z),theta = -35, phi = 30, expand = 0.6, col = "lightgreen",
           ltheta = 120, shade = 0.55, ticktype = "detailed",xlab="Maturity",zlab="Zero-coupon yield",ylab="Time",box=TRUE,border=NA)


