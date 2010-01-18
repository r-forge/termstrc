rm(list = ls())

source("termstrcPackage.R")

## Import CSVs
datadyncouponbonds <- dyncouponbonds(c("DataGermany S.csv", "DataGermany AC.csv", "DataGermany CP.csv"), "GERMANY")

## Nelson/Siegel estimation
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", deltatau = 0.5)
print(ns_res)
summary(ns_res)

## Svensson estimation
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 2)
print(sv_res)
summary(sv_res)







  

