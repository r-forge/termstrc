load("France_1Y.RData")


source("bonddatafunc.R")
source("pricing.R") 
source("rates.R") 

group <- "FRANCE"
dynbonddata  <- dslist

bonddata <- list()
bonddata[[group]] <- dynbonddata[[1]]

data <- prepro_bond("FRANCE",bonddata)

m <-  data$m[[1]]
cf <- data$cf[[1]]

m_p <-  data$m_p[[1]]
cf_p <- data$cf_p[[1]]

beta <- c(0.0511,-0.0124,-0.0303,2.5429)
p <- bond_prices(method="ns",beta,m,cf,lambda)$bond_prices

y <- bond_yields(cf_p,m_p)
