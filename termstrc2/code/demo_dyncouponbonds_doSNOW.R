library("dyntermstrc")

data("GermanBonds")

NODES <- 10 

require("snow")
cl <- makeCluster(NODES, type = "MPI")

require(doSNOW)
registerDoSNOW(cl)

N <- length(datadyncouponbonds)

sv_res <- foreach(i = 1:N, .packages = "dyntermstrc") %dopar% {  
estim_nss(datadyncouponbonds[[i]], c("GERMANY"), method = "sv", deltatau = 2)
}

stopCluster(cl)

class(sv_res) <- "dyntermstrc_nss"

save(sv_res, file = "demo_dyncouponbonds_doSNOW.RData")
