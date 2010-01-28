library("dyntermstrc")

data("GermanBonds")

NODES <- 20 

require("snow")
cl <- makeCluster(NODES, type = "MPI")

require(doSNOW)
registerDoSNOW(cl)

N <- length(datadyncouponbonds)

ns_res <- foreach(i = 1:N, .packages = "dyntermstrc") %dopar% {  
estim_nss(datadyncouponbonds[[i]], c("GERMANY"), method = "ns", deltatau = 0.2)
}

stopCluster(cl)

class(ns_res) <- "dyntermstrc_nss"

save(ns_res, file = "demo_dyncouponbonds_doSNOW3.RData")
