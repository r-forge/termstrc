data(GermanBonds)

## Diebold/Li estimation
dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)

## 3d yield curve plot
plot(dl_res)

## Estimated parameters
plot(param(dl_res))
summary(param(dl_res))

## Estimate Nelson/Siegel model
ns_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", tauconstr = list(c(0.2, 7, 0.2)), optimtype = "allglobal")

## Estimate Svensson model
sv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv",tauconstr = list(c(0.2,7,0.2,0.5)), optimtype = "firstglobal")

## Estimate Adjusted Svensson model
asv_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "asv",tauconstr = list(c(0.2,7,0.2)), optimtype = "firstglobal")
