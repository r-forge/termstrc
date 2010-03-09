oldpar <- par(no.readonly = TRUE)
par(ask=TRUE)

## Load yields matrix from csv
data(zyields)
x <- zyields

maturities <- c(1/12,3/12,6/12,9/12,1:12)
yields <- as.matrix(x[,2:ncol(x)])
dates <- as.Date(x[,1],format="%d.%m.%Y")

## Call class constructor
datazeroyields <- zeroyields(maturities, yields, dates)

## Estimate Diebold/Li model
dl_res <- estim_nss(datazeroyields, "dl", lambda = 1/2)
plot(param(dl_res))

## Estimate Nelson/Siegel model
ns_res <- estim_nss(datazeroyields, "ns")
plot(param(ns_res))

## Estimate Svensson model
sv_res <- estim_nss(datazeroyields, "sv")
plot(param(sv_res))

## Estimate Svensson model with restrictions on the tau parameters
## (this can lead to smoother parameter time series)
sv_res2 <- estim_nss(datazeroyields, "sv", tauconstr =  c(0.2, 3, 0.1,0.5))
plot(param(sv_res2))

## Estimate Adjusted Svensson model
## (this can also lead to smoother parameter time series)
asv_res <- estim_nss(datazeroyields, "asv")
plot(param(asv_res))

## Compare GOF
allgof <- cbind(summary(dl_res)$gof, summary(ns_res)$gof, summary(sv_res)$gof, summary(sv_res2)$gof, summary(asv_res)$gof)
colnames(allgof) <- c("Diebold/Li", "Nelson/Siegel", "Svensson unrestr.", "Svensson", "Adj. Svensson")

par(oldpar)









