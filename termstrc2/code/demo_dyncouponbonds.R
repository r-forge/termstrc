source("termstrcPackage.R")

load("GermanBonds.RData")

## Diebold/Li estimation
dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)

## Nelson/Siegel estimation
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", deltatau = 0.2, optimtype = "allglobal")

## Svensson estimation
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 0.2, optimtype = "allglobal")

## Adjusted Svensson estimation
asv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "asv", deltatau = 0.2, optimtype = "allglobal")

save(datadyncouponbonds, dl_res, ns_res, sv_res, asv_res, file = "demo_dyncoupondbonds_results.RData")



     






  

