rm(list = ls())

load("Germany_040707_030708.RData")
library("rgl")

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


group <- "GERMANY"
method="ns"

matrange <- "all"
myres  <- estim_dyntermstrc(dynbonddata,matrange,method, deltatau = 0.2, diagnosticplots = TRUE)

summary(myres)
myparam <- param(myres)
summary(myparam)
plot(myparam,"3D")
X11()
plot(myparam,"param")
X11()
plot(myparam,"acf") 

method <- "sv"
myres2  <- estim_dyntermstrc(dynbonddata,matrange,method, deltatau = 2, diagnosticplots = TRUE)
