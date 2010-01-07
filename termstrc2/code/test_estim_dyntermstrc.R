load("France_1Y.RData")


source("bonddatafunc.R")
source("create_cf_m.R")
source("pricing.R") 
source("estim_ns.R")
source("helperfunc.R")
source("rates.R")
source("methods_termstrc_ns.R")
source("methods_curves.R")
source("estim_dyntermstrc.R")
source("methods_dyntermstrc.R")
source("param.R")

library(urca) 


group <- "FRANCE"
dynbonddata  <- dslist
method="dl"
fit <- "prices"
weights <- "duration"
matrange <- c(0,20)
myres  <- estim_dyntermstrc(dynbonddata[1:20],matrange,method,fit,weights)

summary(myres)
param <- param(myres)
summary(param)
plot(param,"3D")
plot(param,"diffparam")
plot(param,"acf") 
