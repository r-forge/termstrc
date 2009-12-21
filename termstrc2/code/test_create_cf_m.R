load("France_1Y.RData")


source("bonddatafunc.R")
source("create_cf_m.R")
source("pricing.R") 
source("helperfunc.R")

group <- "FRANCE"
dynbonddata  <- dslist

bonddata <- list()
bonddata[[group]] <- dynbonddata[[1]]

cf <- create_cashflows_matrix(bonddata[[1]])
m  <- create_maturities_matrix(bonddata[[1]])
