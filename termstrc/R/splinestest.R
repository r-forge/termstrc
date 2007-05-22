rm(list = ls())
load("eurobonds.RData")
source("splines.R")
source("tools.R")
source("methods.R")

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- "all"


myres<- splines_estim(group, bonddata, matrange)

#print(myres)
#summary(myres)
#plot(myres)
