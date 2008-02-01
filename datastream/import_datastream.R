
importdatastream <- function(datafile){
#datafile = "austria.csv"
rawdata <- read.csv(datafile, dec=".",sep=",",colClasses = "character")
  data <- list()
  data$ISIN <- rawdata[,"ISIN"]
  data$MATURITYDATE <- as.Date(rawdata[,"MATURITYDATE"],format="%d%m%Y")
  data$ISSUEDATE <- as.Date(rawdata[,"ISSUEDATE"],format="%d%m%Y")
  data$COUPONRATE <- as.numeric(rawdata[,"COUPONRATE"])
  data$PRICE <- as.numeric(rawdata[,"PRICE"])
  data$ACCRUED <- as.numeric(rawdata[,"ACCRUED"])
  data$CASHFLOWS <- list()
  data$TODAY <- as.Date(rawdata[1,"TODAY"],format="%d%m%Y")

  NEXTCOUPON <- ifelse(as.Date(paste(rawdata[,"COUPONDATE"],substring(rawdata[1,"TODAY"],5,8),sep=""),format="%d%m%Y") > data$TODAY,
                       paste(rawdata[,"COUPONDATE"],substring(rawdata[1,"TODAY"],5,8),sep=""),
                       paste(rawdata[,"COUPONDATE"],as.character(as.numeric(substring(rawdata[1,"TODAY"],5,8))+1),sep=""))
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

eurobonds30012008 <- list()
eurobonds30012008$GERMANY <- importdatastream("germany.csv")
eurobonds30012008$AUSTRIA <- importdatastream("austria.csv")

group <- c("GERMANY","AUSTRIA")
bonddata <- eurobonds30012008
matrange <- c(0,20)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "duration"
control <- list(eval.max=100000, iter.max=500)

b <- matrix(c(0,0,0, 1,
			0,0,0, 1),
			nrow=2,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

x <- nelson_estim(group, bonddata, matrange, 
                  method, fit, weights, startparam=b,control)

y <- splines_estim(group, bonddata, matrange)

#print(x)
#summary(x)
#plot(x)
