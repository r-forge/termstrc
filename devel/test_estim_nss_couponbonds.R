rm(list = ls())
library(termstrc)

source("estim_nss_couponbonds.R")
source("spotfwdratedef.R")
source("gradfunc.R")

data(govbonds)

Rprof()
ns_res <- estim_nss(govbonds, c("GERMANY", "AUSTRIA", "FRANCE"), matrange = c(0,30), method = "ns")
Rprof(NULL)
  summaryRprof()


data(GermanBonds)

Rprof()
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv",tauconstr = list(c(0.2,7,0.1,0.5)), optimtype = "firstglobal")
Rprof(NULL)
summaryRprof()

Rprof()
asv_res <- estim_nss(datadyncouponbonds[[1]], c("GERMANY"), method = "asv",tauconstr = list(c(-0.4,10,0.5,0.0)))
Rprof(NULL)
summaryRprof()

## Diebold/Li (German bonds) -> working
dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)

ns_res <- estim_nss(datadyncouponbonds[[1]], c("GERMANY"), method = "ns")

## Nelson/Siegel (German bonds) -> working
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", tauconstr = list(c(1,6,0.1)))
ns_res2 <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", tauconstr = list(c(1,6,0.1)), optimtype = "firstglobal")
matplot(cbind(taustart, param(ns_res)$GERMANY[,4]), type = "l")

pdf("fig_nsbonds.pdf")

taustart <- rep(NA, length(ns_res))
for(i in 1:length(ns_res)) {
  taustart[i] <- ns_res[[i]]$startparam[4]
#plot(ns_res[[i]]$spsearch$GERMANY)
}

dev.off()


