rm(list = ls())
library(termstrc)

source("estim_nss_couponbonds.R")
source("spotfwdratedef.R")
source("gradfunc.R")

data(govbonds)

Rprof()
sv_res2 <- estim_nss(govbonds, c("GERMANY", "AUSTRIA"), matrange = c(0,30), method = "sv")
  Rprof(NULL)
  summaryRprof()


data(GermanBonds)

## Diebold/Li estimation
sv_res <- estim_nss(datadyncouponbonds[[1]], c("GERMANY"), method = "asv",tauconstr = c(0.2,10,0.2,0))
