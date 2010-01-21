rm(list = ls())
source("termstrcPackage.R")
load("demo_dyncoupondbonds_results.RData")
load("GermanBonds.RData")

## Create figures for paper

pdf("fig_sv_dyncouponbonds.pdf", width = 10, height = 7)

op <- par(mfrow = c(2,3))
N <- nrow(param(sv_res2)$GERMANY)
n <- 1

matplot(cbind(param(sv_res)$GERMANY[n:N,1],param(sv_res2)$GERMANY[n:N,1]),type = "l",ylab = expression(hat(beta)[0]),xlab = "Time", lty = c(2,1))
legend("topleft", c("first global", "all global"), lty = c(2,1), col = 1:2)
matplot(cbind(param(sv_res)$GERMANY[n:N,2],param(sv_res2)$GERMANY[n:N,2]),type = "l",ylab = expression(hat(beta)[1]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(sv_res)$GERMANY[n:N,3],param(sv_res2)$GERMANY[n:N,3]),type = "l",ylab = expression(hat(beta)[2]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(sv_res)$GERMANY[n:N,4],param(sv_res2)$GERMANY[n:N,4]),type = "l",ylab = expression(hat(tau)[1]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(sv_res)$GERMANY[n:N,5],param(sv_res2)$GERMANY[n:N,5]),type = "l",ylab = expression(hat(beta)[3]),xlab = "Time", lty = c(2,1))
matplot(cbind(param(sv_res)$GERMANY[n:N,6],param(sv_res2)$GERMANY[n:N,6]),type = "l",ylab = expression(hat(tau)[2]),xlab = "Time", lty = c(2,1))

par(op)
dev.off()

pdf("demo_dyncouponbonds_singlecurves.pdf")

for (i in 1:nrow(param(sv_res2)$GERMANY))
  plot(sv_res2[[i]], matrange = c(0,6))
dev.off()

op <- par(mfrow = c(2,2))
N <- 65#nrow(param(sv_res2)$GERMANY)
n <- 1
for (i in 1:4) {
matplot(cbind(param(ns_res)$GERMANY[n:N,i],param(ns_res2)$GERMANY[n:N,i]),type = "l")
}
par(op)

dev.off()

    plot(param[,1],type="l",xlab="Time",ylab=expression(hat(beta)[0]),
                col=1,lwd=2,... )
           grid()
           plot(param[,2],type="l",xlab="Time",ylab=expression(hat(beta)[1]),
           col=2,lwd=2,... )
           grid()
           plot(param[,3],type="l",xlab="Time",ylab=expression(hat(beta)[2]),
           col=3,lwd=2,... )
           grid()
    
    if(ncol(param)==4) {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
    }
    
    if(ncol(param)==6) {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
           plot(param[,5],type="l",xlab="Time",ylab=expression(hat(beta)[3]),
           col=5,lwd=2,... )
           grid()
           plot(param[,6],type="l",xlab="Time",ylab=expression(hat(tau)[2]),
           col=6,lwd=2,... )
           grid()
    }

