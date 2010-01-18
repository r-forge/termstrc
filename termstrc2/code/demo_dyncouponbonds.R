## Import CSVs
datadyncouponbonds <- dyncouponbonds(c("DataGermany S.csv", "DataGermany AC.csv", "DataGermany CP.csv"), "GERMANY")
datadyncouponbonds <- datadyncouponbonds[10:55]
class(datadyncouponbonds) <- "dyncouponbonds"

## Diebold/Li estimation
dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl")
print(dl_res)
summary(dl_res)

## Nelson/Siegel estimation
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", deltatau = 0.5)
print(ns_res)
plot(ns_res)

tsparam <- param(ns_res)
plot(tsparam)


## Plot startparameter grid search results
plot(ns_res[[1]]$spsearch$GERMANY)

## Svensson estimation
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 2)
print(sv_res)
plot(sv_res)
summary(sv_res)

## Plot startparameter grid search results
plot(sv_res[[1]]$spsearch$GERMANY)








  

