source("termstrcPackage.R")

load("govbonds.RData")
bonddata <- govbonds

test.estim_nss.zeroyields <- function() {
  ## Import CSV
  x <- read.csv("zeroyields.csv",sep=";")

  maturities <- 1:12
  yields <- as.matrix(x[,2:ncol(x)])
  dates <- as.Date(x[,1],format="%d.%m.%Y")

  ## Call class constructor
  datazeroyields <- zeroyields(maturities, yields, dates)
  
  ## Perform Nelson/Siegel estimation
  ns_res <- estim_nss(datazeroyields, "ns", deltatau = 0.2)

  ## Plot startparameters
  plot(ns_res$spsearch)

  ## Plot parameters and curves
  plot(ns_res)
  
  checkTrue(is.numeric(ns_res$optparam))

  ## Perform Svensson estimation
  sv_res <- estim_nss(datazeroyields, "sv", deltatau = 0.2)
  
  ## Plot startparameters
  plot(sv_res$spsearch)

  ## Plot parameters and curves
  plot(sv_res)
  
  checkTrue(is.numeric(sv_res$optparam))
}

test.estim_nss.couponbonds <- function() {
  load("govbonds.RData")
  ns_res <- estim_nss(govbonds, c("GERMANY", "AUSTRIA", "FRANCE"), method = "ns", deltatau = 0.4)
  print(ns_res)
  summary(ns_res)
  checkTrue(is.numeric(ns_res$opt_result$GERMANY$par))
  checkTrue(is.numeric(ns_res$opt_result$AUSTRIA$par))
  checkTrue(is.numeric(ns_res$opt_result$FRANCE$par))
}

test.estim_nss.couponbonds2 <- function() {
  load("govbonds.RData")
  ns_res <- estim_ns(govbonds, c("GERMANY"), method = "ns", deltatau = 1)
  print(ns_res)
  plot(ns_res)
  summary(ns_res)
  plot(param(ns_res))
  dl_res <- estim_ns(govbonds, c("GERMANY"), method = "dl", lambda = 1/ns_res$opt_result$GERMANY$par[4])
  print(dl_res)
  plot(dl_res)
  summary(dl_res)
  plot(param(dl_res))
  
  ## D/L and N/S should have same objective function value
  checkEqualsNumeric(ns_res$opt_result$GERMANY$value, dl_res$opt_result$GERMANY$value)
}

test.estim_nss.couponbonds3 <- function() {
  load("govbonds.RData")
  sv_res <- estim_nss(govbonds, c("GERMANY"), method = "sv", deltatau = 2)
  print(sv_res)
  plot(sv_res)
  plot(param(sv_res))
  summary(sv_res)
  checkTrue(is.numeric(sv_res$opt_result$GERMANY$par))
}






