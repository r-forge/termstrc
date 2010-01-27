rm(list = ls())
source("termstrcPackage.R")
load("demo_dyncoupondbonds_asv_results.RData")
load("GermanBonds.RData")

## Create figures for paper

pdf("fig_asv_dyncouponbonds.pdf", width = 10, height = 7)

op <- par(mfrow = c(2,3))
N <- nrow(param(asv_res2)$GERMANY)
n <- 1

matplot(cbind(param(asv_res)$GERMANY[n:N,1],param(asv_res2)$GERMANY[n:N,1]),type = "l",ylab = expression(hat(beta)[0]),xlab = "Time", lty = c(2,1))
legend("topleft", c("first global", "all global"), lty = c(2,1), col = 1:2)
matplot(cbind(param(asv_res)$GERMANY[n:N,2],param(asv_res2)$GERMANY[n:N,2]),type = "l",ylab = expression(hat(beta)[1]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(asv_res)$GERMANY[n:N,3],param(asv_res2)$GERMANY[n:N,3]),type = "l",ylab = expression(hat(beta)[2]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(asv_res)$GERMANY[n:N,4],param(asv_res2)$GERMANY[n:N,4]),type = "l",ylab = expression(hat(tau)[1]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(asv_res)$GERMANY[n:N,5],param(asv_res2)$GERMANY[n:N,5]),type = "l",ylab = expression(hat(beta)[3]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(asv_res)$GERMANY[n:N,6],param(asv_res2)$GERMANY[n:N,6]),type = "l",ylab = expression(hat(tau)[2]),xlab = "Time", lty = c(2,1))

par(op)
dev.off()
