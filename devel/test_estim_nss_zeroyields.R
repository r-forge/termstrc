rm(list = ls())
library("termstrc")

source("spotfwdratedef.R")
source("estim_nss_zeroyields.R")

data(zyields)
x <- zyields

maturities <- c(1/12,3/12,6/12,9/12,1:12)
yields <- as.matrix(x[,2:ncol(x)])
dates <- as.Date(x[,1],format="%d.%m.%Y")

## Call class constructor

datazeroyields <- zeroyields(maturities, yields, dates)

dl_res <- estim_nss(datazeroyields, "dl")
ns_res <- estim_nss(datazeroyields, "ns")

sv_res <- estim_nss(datazeroyields, "sv", tauconstr =  c(0.2, 3, 0.1,0.5))
sv_res2 <- estim_nss(datazeroyields, "sv")

asv_res <- estim_nss(datazeroyields, "asv", tauconstr =  c(0.1, 4, 0.1))
asv_res2 <- estim_nss(datazeroyields, "asv")

# Compare GOF

allgof <- cbind(summary(dl_res)$gof, summary(ns_res)$gof, summary(sv_res)$gof, summary(sv_res2)$gof, summary(asv_res)$gof, summary(asv_res2)$gof)
colnames(allgof) <- c("Diebold/Li", "Nelson/Siegel", "Svensson", "Svensson unrestr.", "Adj. Svensson", "Adj. Svensson unrestr.")

pdf("allparam.pdf")
plot(param(dl_res))
plot(param(ns_res))
plot(param(sv_res))
plot(param(sv_res2))
plot(param(asv_res))
plot(param(asv_res2))
dev.off()

pdf("fcontribasv.pdf")
for (i in 1:80) fcontrib(param(asv_res), method = "asv", index = i)
dev.off()
pdf("fcontribasv.pdf")
for (i in 1:80) fcontrib(param(asv_res2), method = "asv", index = i)
dev.off()



