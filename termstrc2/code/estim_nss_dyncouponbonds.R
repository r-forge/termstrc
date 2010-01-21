#############################################################################
### Nelson/Siegel-type yield curve estimation method for 'dyncouponbonds' ###
#############################################################################

estim_nss.dyncouponbonds <- function(dynbonddata, group, matrange="all",method="ns",
                              lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                              deltatau=1,              # interval for parameter grid
                              optimtype = "firstglobal", # 'firstglobal' of 'allglobal'
                             # nlmbinbOptions = list(control = list(iter.max = 1500, eval.max = 1500,trace)),
                              constrOptimOptions = list(control = list(maxit = 2000), outer.iterations = 200, outer.eps = 1e-04)
                     ) {
  res <- list()

  ## perform sequence of term structure estimations
  for (i in seq(length(dynbonddata))) {
    if(i>1 && optimtype == "firstglobal"){
      ## use optimal parameters from previous period as start parameters
      b <- t(mapply(function(j) res[[i-1]]$opt_result[[j]]$par,  seq_along(group)))
      rownames(b) <- group                               
    } else b <- NULL
    
    ## static estimation
    bonddata <- dynbonddata[[i]]
    res[[i]] <- estim_nss(bonddata, group, matrange, method=method, startparam=b, lambda=lambda,deltatau,constrOptimOptions)
  }
  class(res) <- "dyntermstrc_nss"

  res
}
