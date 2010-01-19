load("govbonds.RData")
ns_res <- estim_nss(govbonds, c("GERMANY"),matrange = c(0,30), method = "ns", deltatau = 2)
ns_res <- estim_nss(govbonds, c("GERMANY", "AUSTRIA", "FRANCE"),matrange = c(0,30), method = "ns", deltatau = 1)
print(ns_res)
plot(ns_res)
summary(ns_res)

## Plot startparameter grid search results
par(mfrow=c(1,3))
plot(ns_res$spsearch$GERMANY,main="GERMANY")
plot(ns_res$spsearch$AUSTRIA,main="AUSTRIA")
plot(ns_res$spsearch$FRANCE,main="FRANCE")

## Plot all yield curves in one figure  
par()
plot(ns_res,multiple=TRUE)


## Svensson estimation with parameter grid search (can take some time, increase deltatau for wider grid)
sv_res <- estim_nss(govbonds, c("GERMANY"),matrange = c(0,30), method = "sv", deltatau = 1)
print(sv_res)
plot(sv_res)
summary(sv_res)

## Plot startparameter grid search results
plot(sv_res$spsearch$GERMANY)

## Cubic splines estimation
cs_res <- estim_cs(govbonds,c("GERMANY"),matrange=c(0,30))
print(cs_res)
summary(cs_res)
plot(cs_res)
## Pricing errors per bond
plot(cs_res,errors="price",inset=c(.1,.3))

## Compare model performance (gofs)
