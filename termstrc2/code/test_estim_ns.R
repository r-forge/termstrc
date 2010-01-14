rm(list = ls())

library("rgl")
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

dl_res <- estim_ns(govbonds, c("GERMANY"), method = "dl", lamda = 1/5)

ns_res <- estim_ns(govbonds, c("GERMANY"), "ns", deltatau = 0.5)


#sv_res <- estim_ns(govbonds, c("GERMANY"), method = "sv", deltatau = 5)

