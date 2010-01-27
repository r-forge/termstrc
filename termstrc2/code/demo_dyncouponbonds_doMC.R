source("termstrcPackage.R")
load("GermanBonds.RData")

library(doMC)
registerDoMC()
workers <- getDoParWorkers()

N <- length(datadyncouponbonds)
N <- 4

sv_res <- foreach(i = 1:N) %dopar% { 
  sv_res <- estim_nss(datadyncouponbonds[[i]], c("GERMANY"), method = "sv", deltatau = 0.2)
}

class(sv_res) <- "dyntermstrc_nss"

save(workers, sv_res, file = "demo_dyncouponbonds_doMC.RData")
