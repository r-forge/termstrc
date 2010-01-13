
## dynamic estimation of the term structure

estim_dyntermstrc <- function(dynbonddata,matrange="all",method="ns",
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
      b <- switch(method,
                  "ns" = matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=4,byrow=TRUE),
                  "sv" =matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=6,byrow=TRUE),
                  "dl" =matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=3,byrow=TRUE))
      rownames(b) <- group
      colnames(b) <- switch(method,
                            "ns" = c("beta0","beta1","beta2","tau1"),
                            "sv" = c("beta0","beta1","beta2","tau1","beta3","tau2"),
                            "dl" = c("beta0","beta2","beta3"))
                            
    } else b <- NULL

    # static estimation
    group <- names(dynbonddata)[i]
    bonddata <- list()
    bonddata[[group]] <- dynbonddata[[i]]
    res[[i]] <- estim_ns(bonddata=bonddata,group, matrange, method=method, startparam=b, lambda=lambda,deltatau,control,outer.iterations,outer.eps)
  }
  class(res) <- "dyntermstrc"

  res
}
