load("France_1Y.RData")


source("bonddatafunc.R")
source("create_cf_m.R")
source("pricing.R") 
source("helperfunc.R")

group <- "FRANCE"
dynbonddata  <- dslist

bonddata <- list()
bonddata[[group]] <- dynbonddata[[1]]



# remove bonds with ISIN "FR0000571218","FR0108354806"  from bonddata
newbonddata <- rm_bond(bonddata,c("FR0000571218","FR0108354806"),group)

# remove bonds with ISIN "FR0000571218","FR0108354806"  from dynbonddata
newdynbonddata <- dyn_rm_bond(dynbonddata,c("FR0000571218","FR0108354806"))

# preprocess bonddata, i.e., calculate maturity, cashflowmatrix, sort data, calculate yield etc. 
data <- prepro_bond("FRANCE",bonddata)

data$m[[1]]
data$cf[[1]]



