rm(list = ls())

source("methods.R")
source("nelson.R")
source("splines.R")
source("tools.R")

load("eurobonds.RData")

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- c(2,10) 

x <- splines_estim(group, bonddata, matrange)

print(x)
summary(x)
plot(x)
