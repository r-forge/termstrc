source("termstrcPackage.R")

load("GermanBonds.RData")

## Svensson estimation
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 0.5, optimtype = "allglobal")
sv_res2 <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 0.5, optimtype = "firstglobal")

dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)

ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", deltatau = 0.2, optimtype = "allglobal")
ns_res2 <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", deltatau = 0.2, optimtype = "allglobal")

save(sv_res, sv_res2, dl_res, ns_res, ns_res2, file = "demo_dyncoupondbonds_results.RData")



     






  

