rm(list = ls())

load("Germany_040707_030708.RData")

source("bonddatafunc.R")
source("create_cf_m.R")
source("pricing.R") 
source("cubicfunc.R")
source("estim_ns.R")
source("helperfunc.R")
source("rates.R")
source("methods_termstrc_ns.R")
source("methods_curves.R")

sv_res <- estim_ns(dynbonddata, c("GERMANY"), method = "sv", deltatau = 0.2, diagnosticplots = FALSE)
save(sv_res, file = "sv_res.RData")
