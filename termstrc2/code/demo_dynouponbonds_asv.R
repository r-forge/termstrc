source("termstrcPackage.R")

load("GermanBonds.RData")

## Adjusted Svensson estimation
asv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "asv", deltatau = 0.5, optimtype = "allglobal")
asv_res2 <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "asv", deltatau = 0.5, optimtype = "firstglobal")

save(asv_res, asv_res2, file = "demo_dyncoupondbonds_asv_results.RData")
