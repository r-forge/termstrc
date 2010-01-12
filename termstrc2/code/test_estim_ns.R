rm(list = ls())

load("govbonds.RData")

source("bonddatafunc.R")
source("create_cf_m.R")
source("pricing.R") 
source("cubicfunc.R")
source("estim_ns.R")
source("helperfunc.R")
source("rates.R")
source("methods_termstrc_ns.R")
source("methods_curves.R")

ns_res <- estim_ns(govbonds, c("GERMANY"), "ns", deltatau = 0.4, diagnosticplots = TRUE)
