data(GermanBonds)

## Diebold/Li estimation
dl_res <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "dl", lambda = 1/3)

## 3d yield curve plot
plot(dl_res)

## Estimated parameters
plot(param(dl_res))
summary(param(dl_res))
