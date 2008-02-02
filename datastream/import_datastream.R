
importdatastream <- function(datafile){
#datafile = "belgium.csv"
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

cf_p <- create_cashflows_matrix(data,include_price=T)
m_p <- create_maturities_matrix(data,include_price=T)


govbonds <- list()
govbonds$GERMANY <- importdatastream("germany.csv")
govbonds$AUSTRIA <- importdatastream("austria.csv")
govbonds$BELGIUM <- importdatastream("belgium.csv")
govbonds$FINLAND <- importdatastream("finland.csv")
govbonds$FRANCE <- importdatastream("france.csv")
govbonds$SPAIN <- importdatastream("spain.csv")

group <- c("GERMANY","AUSTRIA","BELGIUM","FINLAND","FRANCE","SPAIN")

bonddata <- govbonds
matrange <- c(1,20)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "duration"
control <- list(eval.max=100000, iter.max=500)

b <- matrix(c(0.1,0.1,0.1, 1,
	      0.1,0.1,0.1, 1,
              0.1,0.1,0.1, 1,
              0.1,0.1,0.1, 1,
              0.1,0.1,0.1, 1,
	      0.1,0.1,0.1, 1),
	      nrow=6,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

x <- nelson_estim(group, bonddata, matrange, 
                  method, fit, weights, startparam=b,control)

group <- c("GERMANY","AUSTRIA","SPAIN")
y <- splines_estim(group, bonddata, matrange)

#print(x)
#summary(x)
#plot(x)
