print.dyntermstrc_yields <- function(x, ...){
  cat("---------------------------------------------------\n")
  cat("Parameters for yield based dynamic term structure estimation:\n")
  cat("---------------------------------------------------\n")
  cat("Method:",get_realnames(x$method),"\n")
  cat("Number of oberservations:",nrow(x$optparam),"\n")
  cat("---------------------------------------------------\n")
  cat("Parameter summary:\n")
  cat("---------------------------------------------------\n")
  tsparam <- param.dyntermstrc_yields(x)
  print(lapply(tsparam,summary.default))
  cat("---------------------------------------------------\n")
}



summary.dyntermstrc_yields <- function(object, ...){
  x <- object
  y_mrsme <-  round(mean(sqrt(apply((x$yields-x$yhat)^2,2,mean))),6)
  y_maabse <- round(mean(apply(abs(x$yields-x$yhat),2,mean)),6)
  sumry <- list()
  sumry$gof <- rbind(y_mrsme,y_maabse)
  rownames(sumry$gof) <- c("mean RMSE-Yields (in %)", "mean AABSE-Yields (in %)")

  if (object$method != "dl") {
    ## extract convergence info
    for (i in (1:length(x$opt_result))) {
      sumry$convergence[i] <- x$opt_result[[i]]$convergence
      sumry$solvermsg <- x$opt_result[[i]]$message
    }
  }  
  class(sumry) <- "summary.dyntermstrc_yields"
  sumry
}


print.summary.dyntermstrc_yields <- function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")

    print.default(x$gof)

    if (length(x) > 1) {
      cat("\n")
      cat("---------------------------------------------------\n")
      cat("Convergence information from optim ():\n")
      cat("---------------------------------------------------\n")
      
      print.default(x$convergence)
    }
  }

plot.dyntermstrc_yields <- function(x,...)
  {
    ## plot estimated yield curves in 3D
 
    Z <- matrix(nrow=nrow(x$optparam),ncol=length(x$maturities))# OK
    for (i in 1:nrow(x$optparam)){
      Z[i,] <- spotrates(x$method,x$optparam[i,],x$maturities, x$lambda)
    }

    X <- 1:nrow(Z)
    Y <- x$maturities
    
    open3d()
    persp3d(X, Y, Z, col = "green3", box = FALSE,xlab = "Dates", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
    
  }


