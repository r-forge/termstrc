rm(list = ls())
load("eurobonds.RData")
source("splines.R")
source("tools.R")
source("methods.R")

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
maturity_spectrum <- "all"


myres<- splines_estim(group, bonddata, maturity_spectrum)

#print(myres)
#summary(myres)
#plot(myres)
