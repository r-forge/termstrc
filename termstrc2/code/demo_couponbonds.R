load("govbonds.RData")
ns_res <- estim_nss(govbonds, c("GERMANY", "AUSTRIA", "FRANCE"),matrange = c(0,30), method = "ns", deltatau = 0.2)
print(ns_res)
plot(ns_res)
summary(ns_res)

## Svensson estimation with parameter grid search (can take some time, increase deltatau for wider grid)
sv_res <- estim_nss(govbonds, c("GERMANY"),matrange = c(0,30), method = "sv", deltatau = 1)
print(sv_res)
plot(sv_res)
summary(sv_res)

## Plot startparameter grid search results
plot(sv_res$spsearch$GERMANY)

