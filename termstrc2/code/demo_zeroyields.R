rm(list = ls())

source("termstrcPackage.R")
#mtrace(estim_nss.zeroyields)
#mtrace(estimateyieldcurve)

## Import CSV

x <- read.csv("zeroyields.csv",sep=";")

maturities <- 1:12
yields <- as.matrix(x[100:150,2:13])
dates <- as.Date(x[100:150,1],format="%d.%m.%Y")

## Call class constructor

datazeroyields <- zeroyields(maturities, yields, dates)

## Perform Nelson/Siegel estimation
#ns_res <- estim_nss(datazeroyields, "ns", deltatau = 0.2)

## Plot startparameters
#plot(ns_res$spsearch)

## Perform Svensson estimation

sv_res <- estim_nss(datazeroyields, "sv", deltatau = 0.2)


## Plot startparameters
plot(sv_res$spsearch)





