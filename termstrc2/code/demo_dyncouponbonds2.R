source("termstrcPackage.R")

datadyncouponbonds <- dyncouponbonds(c("DataGermany S.csv", "DataGermany AC.csv", "DataGermany CP.csv"), "GERMANY")

dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", deltatau = 0.2, optimtype = "allglobal")
save(dl_res, ns_res, file = "demo_dyncoupondbonds2_results.RData")


     






  

