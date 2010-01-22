rm(list = ls())
source("termstrcPackage.R")
## Import CSV
x <- read.csv("zeroyields.csv",sep=";")

maturities <- c(1/12,3/12,6/12,9/12,1:12)
yields <- as.matrix(x[,2:ncol(x)])
dates <- as.Date(x[,1],format="%d.%m.%Y")

## Call class constructor

datazeroyields <- zeroyields(maturities, yields, dates)

## Perform Nelson/Siegel estimation
ns_res <- estim_nss(datazeroyields, "ns", deltatau = 0.2, optimtype = "firstglobal")
ns_res2 <- estim_nss(datazeroyields, "ns", deltatau = 0.2, optimtype = "allglobal")

## Plot startparameters
plot(ns_res$spsearch)

## Plot parameters and curves
plot(ns_res)

## Perform Svensson estimation
sv_res <- estim_nss(datazeroyields, "sv", deltatau = 0.2)
sv_res2 <- estim_nss(datazeroyields, "sv", deltatau = 0.2, optimtype = "allglobal")

## Perform Ajusted Svensson estimation
asv_res <- estim_nss(datazeroyields, "asv", deltatau = 0.2)
asv_res2 <- estim_nss(datazeroyields, "asv", deltatau = 0.2, optimtype = "allglobal")

## Plot startparameters
plot(sv_res$spsearch)

## Plot parameters and curves
plot(sv_res)

## Plot single curves
pdf("zerocurves.pdf")
yhat_ns= matrix(nrow=nrow(ns_res$optparam),ncol=length(ns_res$maturities))
yhat_sv= matrix(nrow=nrow(sv_res$optparam),ncol=length(sv_res$maturities))
for (i in 1:nrow(ns_res$yields)) {
  plot(ns_res$maturities, ns_res$yields[i,])
  yhat_ns[i,] <- spr_ns(ns_res$optparam[i,],ns_res$maturities)
  points(ns_res$maturities, yhat_ns[i,], col = "red")
  yhat_sv[i,] <- spr_sv(sv_res$optparam[i,],sv_res$maturities)
  points(sv_res$maturities, yhat_sv[i,], col = "blue") 
}
 
dev.off()





