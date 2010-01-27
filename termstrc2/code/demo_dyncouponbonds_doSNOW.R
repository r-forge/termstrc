source("termstrcPackage.R")
load("GermanBonds.RData")

NODES <- 5 

library("snow")
cl <- makeCluster(NODES, type = "MPI")

library(doSNOW)
registerDoSNOW(cl)

N <- length(datadyncouponbonds)
N <- 10

sv_res <- foreach(i = 1:N) %dopar% { 
  sv_res <- estim_nss(datadyncouponbonds[[i]], c("GERMANY"), method = "sv", deltatau = 1)
}

stopCluster(cl)

class(sv_res) <- "dyntermstrc_nss"

save(sv_res, file = "demo_dyncouponbonds_doSNOW.RData")
