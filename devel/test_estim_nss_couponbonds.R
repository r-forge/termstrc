rm(list = ls())
library(termstrc)

source("estim_nss_couponbonds.R")
source("spotfwdratedef.R")
source("gradfunc.R")

data(govbonds)

Rprof()
sv_res2 <- estim_nss(govbonds, c("GERMANY"), matrange = c(0,30), method = "sv")
  Rprof(NULL)
  summaryRprof()


data(GermanBonds)

sv_res <- estim_nss(datadyncouponbonds[[1]], c("GERMANY"), method = "sv",tauconstr = c(1,10,0.2,0.5))

seq(tauconstr[1] + tauconstr[3], tauconstr[2] - tauconstr[3], tauconstr[3])
