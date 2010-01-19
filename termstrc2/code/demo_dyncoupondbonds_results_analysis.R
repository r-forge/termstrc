rm(list = ls())
source("termstrcPackage.R")
load("demo_dyncoupondbonds_results.RData")
load("demo_dyncoupondbonds2_results.RData")

## Import CSVs
datadyncouponbonds <- dyncouponbonds(c("DataGermany S.csv", "DataGermany AC.csv", "DataGermany CP.csv"), "GERMANY")

## Nelson/Siegel estimation
ns_res2 <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "ns", deltatau = 1, optimtype = "firstglobal")

## Svensson estimation
sv_res2 <- estim_nss(datadyncouponbonds, c("GERMANY"), method = "sv", deltatau = 1, optimtype = "firstglobal")

summary(sv_res)
summary(sv_res2)

param_names <- c("beta_0", "beta_1","beta_2","tau_1","beta_3","tau_2")

pdf("demo_dyncouponbonds_singlecurves.pdf")

for (i in 1:length(param(sv_res2)$GERMANY[,i]))
  plot(sv_res2[[i]])
dev.off()

pdf("demo_dyncouponbonds_results.pdf")
#plot(param(dl_res))

op <- par(mfrow = c(2,3))
N <- 55#length(param(sv_res)$GERMANY[,i])
n <- 8
for (i in 1:6) {
matplot(cbind(param(sv_res)$GERMANY[n:N,1],param(sv_res2)$GERMANY[n:N,1]),type = "l",ylab = expression(hat(beta)[0]),xlab = "Time")
matplot(cbind(param(sv_res)$GERMANY[n:N,2],param(sv_res2)$GERMANY[n:N,2]),type = "l",ylab = expression(hat(beta)[1]),xlab = "Time")
matplot(cbind(param(sv_res)$GERMANY[n:N,3],param(sv_res2)$GERMANY[n:N,3]),type = "l",ylab = expression(hat(beta)[2]),xlab = "Time")
matplot(cbind(param(sv_res)$GERMANY[n:N,4],param(sv_res2)$GERMANY[n:N,4]),type = "l",ylab = expression(hat(tau)[1]),xlab = "Time")
matplot(cbind(param(sv_res)$GERMANY[n:N,5],param(sv_res2)$GERMANY[n:N,5]),type = "l",ylab = expression(hat(beta)[3]),xlab = "Time")
matplot(cbind(param(sv_res)$GERMANY[n:N,6],param(sv_res2)$GERMANY[n:N,6]),type = "l",ylab = expression(hat(tau2)[2]),xlab = "Time")
}
par(op)

op <- par(mfrow = c(2,2))
N <- 55#length(param(ns_res)$GERMANY[,i])
n <- 8
for (i in 1:4) {
matplot(cbind(param(ns_res)$GERMANY[n:N,i],param(ns_res2)$GERMANY[n:N,i]),type = "l",ylab = param_names[i])
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

