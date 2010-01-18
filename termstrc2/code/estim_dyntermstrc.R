
## dynamic estimation of the term structure

estim_nss.dyncouponbonds <- function(dynbonddata, group, matrange="all",method="ns",
                              lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                              deltatau=1,              # interval for parameter grid
                              control=list(),            # options or optim() 
                              outer.iterations = 30,     # options for constrOptim()
                              outer.eps = 1e-04
                     ) {
  res <- list()

 
  # perform sequence of term structure estimations
  for (i in seq(length(dynbonddata))) {
    if(i>1){
      # use optimal parameters from previous period as start parameters
      b <- t(mapply(function(j) res[[i-1]]$opt_result[[j]]$par,  seq_along(group)))
      rownames(b) <- group                               
    } else b <- NULL

    # static estimation
    res[[i]] <- estim_nss(bonddata=dynbonddata[[i]],group, matrange, method=method, startparam=b, lambda=lambda,deltatau,control,outer.iterations,outer.eps)
  }
  class(res) <- "dyntermstrc_nss"

  res
}
