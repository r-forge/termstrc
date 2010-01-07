
## dynamic estimation of the term structure

estim_dyntermstrc <- function(dynbonddata,matrange,method,fit,weights,
                        startparam="auto",lambda=0.0609,otype="nlminb",...) {

res <- list()
 
 
  # perform sequence of term structure estimations
  for (i in seq(length(dynbonddata))) {
    if(i>1){
      # use optimal parameters from previous period as start parameters
      b <- switch(method,
                  "ns" = matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=4,byrow=TRUE),
                  "sv" =matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=6,byrow=TRUE),
                  "dl" =matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=3,byrow=TRUE))
      rownames(b) <- group
      colnames(b) <- switch(method,
                            "ns" = c("beta0","beta1","beta2","tau1"),
                            "sv" = c("beta0","beta1","beta2","tau1","beta3","tau2"),
                            "dl" = c("beta0","beta2","beta3"))
                            
    } else b <- "auto"

    # static estimation
    group <- names(dynbonddata)[i]
    bonddata <- list()
    bonddata[[group]] <- dynbonddata[[i]]
    res[[i]] <- estim_ns(bonddata=bonddata,group, matrange, 
                           method=method, fit, weights, startparam=b,
                           lambda=lambda,otype=otype,...)
  }
  class(res) <- "dyntermstrc"

  res
}
