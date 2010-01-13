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
myres  <- estim_dyntermstrc(dynbonddata,matrange,method, deltatau = 0.2, diagnosticplots = FALSE)

method <- "sv"
myres2  <- estim_dyntermstrc(dynbonddata,matrange,method, deltatau = , diagnosticplots = FALSE)

save.image(file = "test_cluster_estim_dyntermstrc.RData")
