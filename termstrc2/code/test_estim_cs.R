load("France_1Y.RData")


source("bonddatafunc.R")
source("create_cf_m.R")
source("pricing.R") 
source("cubicfunc.R")
source("estim_cs.R")
source("helperfunc.R")
source("rates.R")
source("methods_termstrc_cs.R")
source("methods_curves.R")

library(robustbase)
library(sandwich)

group <- "FRANCE"
dynbonddata  <- dslist

bonddata <- list()
bonddata[[group]] <- dynbonddata[[1]]


cs_res <- estim_cs(bonddata,group,rse=FALSE)
# adapt summary method for ks.test 
