summary.dyntermstrc <- function(object, ...) {
  x <- object

  # extract convergence info
  sumry <- list()
  sumry$convergence <- vector()
  sumry$solvermsg <- vector()
  for (i in 1:length(x)) {
    sumry$convergence[i] <- summary(x[[i]])$convergencegroup[[1]]
    sumry$solvermsg[i] <- summary(x[[i]])$convergence[[1]]
  }

  sumry$convprobs <- which(sumry$converge != "converged")
  class(sumry) <- "summary.dyntermstrc"
  sumry
}



print.summary.dyntermstrc <- function(x,...) {
    cat("Convergence info:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    if(length(x$convprobs)==0) {
      print.default("No convergence problems!")
    } else {
      cat("Convergence problems are at:\n")
      print.default(x$convprobs)
    }
    cat("---------------------------------------------------\n")
}




plot.dyntermstrc <- function(x,range=c(0,20), ...) {

  old.par <- par(no.readonly = TRUE) 
  par(mfrow=c(1,1),...)
      
  # 3D plot of zero-coupon yield curves
    param <- param.dyntermstrc(x)

    X <- seq(if(range[1]==0) range[1]+0.1 else range[1],range[2],0.1)
    Y <- seq(nrow(param))

    Z <-  mapply(function(i) spotrates(method=x[[1]]$method,param[i,],X,x[[1]]$lambda), seq(nrow(param)))
    open3d()
    persp3d(X,Y,Z,col = "green3",xlab="Maturity (years)",
          zlab="Zero-coupon yields (in %)",ylab="Time",box=FALSE)

 
on.exit(par(old.par))
}
