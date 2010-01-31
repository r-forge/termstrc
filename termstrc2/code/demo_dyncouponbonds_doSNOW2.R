library("dyntermstrc")

data("GermanBonds")

NODES <- 20 

require("snow")
cl <- makeCluster(NODES, type = "MPI")

require(doSNOW)
registerDoSNOW(cl)

N <- length(datadyncouponbonds)

asv_res <- foreach(i = 1:N, .packages = "dyntermstrc") %dopar% {  
estim_nss(datadyncouponbonds[[i]], c("GERMANY"), method = "asv", deltatau = 0.1)
}

stopCluster(cl)

class(asv_res) <- "dyntermstrc_nss"

save(asv_res, file = "demo_dyncouponbonds_doSNOW2.RData")
