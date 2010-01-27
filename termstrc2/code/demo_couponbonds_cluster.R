rm(list = ls())
source("termstrcPackage.R")
load("GermanBonds.RData")

sv_res <- estim_nss(datadyncouponbonds[[1]], c("GERMANY"), method = "sv", deltatau = 0.2, optimtype = "allglobal")

save(sv_res, file = "demo_dyncouponbonds_cluster.RData")
