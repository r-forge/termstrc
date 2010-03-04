rm(list = ls())
library(termstrc)

source("estim_nss_couponbonds.R")
source("spotfwdratedef.R")
source("gradfunc.R")

data(govbonds)

Rprof()
sv_res2 <- estim_nss(govbonds, c("GERMANY"), matrange = c(0,30), method = "dl")
Rprof(NULL)
  summaryRprof()


data(GermanBonds)

Rprof()
sv_res <- estim_nss(datadyncouponbonds[[1]], c("GERMANY"), method = "sv",tauconstr = list(c(-0.4,10,0.5,0.1)))
Rprof(NULL)
summaryRprof()

# Diebold/Li (German bonds) -> working
dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)

# Nelson/Siegel (German bonds) -> working
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", tauconstr = list(c(1,6,0.2,0)))



