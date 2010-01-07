load("France_1Y.RData")


source("bonddatafunc.R")
source("create_cf_m.R")
source("pricing.R") 
source("cubicfunc.R")
source("estim_ns.R")
source("helperfunc.R")
source("rates.R")
source("methods_termstrc_ns.R")
source("methods_curves.R")


group <- "FRANCE"
dynbonddata  <- dslist
method="ns"

bonddata <- list()
bonddata[[group]] <- dynbonddata[[1]]


cs_res <- estim_ns(bonddata,group,method=method)
