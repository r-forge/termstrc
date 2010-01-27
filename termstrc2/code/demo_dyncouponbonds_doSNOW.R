source("termstrcPackage.R")
load("GermanBonds.RData")

NODES <- 5 

require("snow")
cl <- makeCluster(NODES, type = "MPI")

require(doSNOW)
registerDoSNOW(cl)

N <- length(datadyncouponbonds)
N <- 5

sv_res <- foreach(i = 1:N, packages = c("rgl", "urca", "sandwich")) %dopar% { 
  ## Estimation kernel
source("estim_cs.R")
source("estim_nss_couponbonds.R")
source("estim_nss_dyncouponbonds.R")
source("estim_nss_zeroyields.R")
source("bondpricing.R")

## Parametric forms of spot rate functions
source("spotfwdratedef.R")
source("cubicfunc.R")

## Data handling
source("couponbonds_data.R")
source("create_cf_m.R")
source("methods_couponbonds.R")

## Methods for estimation results
source("methods_curves.R")
source("methods_dyntermstrc_nss.R")
source("methods_dyntermstrc_param.R")
source("methods_dyntermstrc_yields.R")
source("methods_termstrc_cs.R")
source("methods_termstrc_nss.R")
source("methods_zeroyields.R")
source("gof.R")
source("factorcontrib.R")  
estim_nss(datadyncouponbonds[[i]], c("GERMANY"), method = "sv", deltatau = 1)
}

stopCluster(cl)

class(sv_res) <- "dyntermstrc_nss"

schas = 2,

save(schas,sv_res, file = "demo_dyncouponbonds_doSNOW.RData")
