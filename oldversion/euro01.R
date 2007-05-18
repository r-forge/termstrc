#data(eurobonds)

rm(list=ls())
source("svensson.R")
source("tools.R")
load("eurobonds.RData")

countries <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
maturity_spectrum <- "all"
method = "Nelson/Siegel"
fit = "prices"
weights = "none"
control=list(eval.max=100000)

b<-matrix(c(0.02547394, -0.012162592, -0.02547394,    1,
 			0.02611532, -0.011367422, -0.02611532,    1,
			0.02578871, -0.015207250, -0.02578871,    1),
			nrow=3,ncol=4,byrow=TRUE)
			
rownames(b)<-countries

colnames(b)<-c("beta0","beta1","beta2","tau1")

myres<- termstrc_estim(countries, bonddata, maturity_spectrum, 
    method, fit, weights, startparam=b,control)



    
#print(myres)
#summary(myres)
#plot(myres)