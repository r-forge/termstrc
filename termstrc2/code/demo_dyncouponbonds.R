## Import CSVs
datadyncouponbonds <- dyncouponbonds(c("DataGermany S.csv", "DataGermany AC.csv", "DataGermany CP.csv"), "GERMANY")

## Svensson estimation
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 2)
print(sv_res)
plot(sv_res)
summary(sv_res)

## Plot startparameter grid search results
plot(sv_res[[1]]$spsearch$GERMANY)

## Plot parameters

plot(param(sv_res))
summary(param(sv_res))

     






  

