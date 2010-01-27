data(zeroyields)
x <- zeroyields

maturities <- c(1/12,3/12,6/12,9/12,1:12)
yields <- as.matrix(x[,2:ncol(x)])
dates <- as.Date(x[,1],format="%d.%m.%Y")

## Call class constructor

datazeroyields <- zeroyields(maturities, yields, dates)

## Perform Diebold/Li estimation
dl_res <- estim_nss(datazeroyields, "dl", lambda = 1/2)
summary(dl_res)

## Perform Nelson/Siegel estimation
ns_res <- estim_nss(datazeroyields, "ns", deltatau = 0.2, optimtype = "allglobal")
summary(ns_res)

## Plot startparameters
plot(ns_res$spsearch[[1]])

## Plot parameters and curves
plot(ns_res)

## Perform Svensson estimation
sv_res <- estim_nss(datazeroyields, "sv", deltatau = 0.2)
sv_res2 <- estim_nss(datazeroyields, "sv", deltatau = 0.2, optimtype = "allglobal")

## Perform Adjusted Svensson estimation
asv_res <- estim_nss(datazeroyields, "asv", deltatau = 0.2)
asv_res2 <- estim_nss(datazeroyields, "asv", deltatau = 0.2, optimtype = "allglobal")

## Plot startparameters
plot(sv_res$spsearch[[1]])

## Plot parameters and curves
plot(sv_res)






