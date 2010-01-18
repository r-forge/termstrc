source("termstrcPackage.R")

## Import CSVs
datadyncouponbonds <- dyncouponbonds(c("DataGermany S.csv", "DataGermany AC.csv", "DataGermany CP.csv"), "GERMANY")

## datadyncouponbonds <- datadyncouponbonds[1:3]
## class(datadyncouponbonds) <- "dyncouponbonds"

## Svensson estimation
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 1, optimtype = "allglobal")

save(sv_res, file = "demo_dyncoupondbonds_results.RData")

## print(sv_res)
## plot(sv_res)
## summary(sv_res)

## ## Plot startparameter grid search results
## plot(sv_res[[1]]$spsearch$GERMANY)

## ## Plot parameters

## plot(param(sv_res))
## summary(param(sv_res))

     






  

