### Nelson/Siegel and Svensson spot curve estimaton from zero yields 

zeroyields <- function(maturities, yields, dates)
  {
    zy <- list(maturities = maturities, yields = yields, dates = dates)
    class(zy) <- "zeroyields"
    zy
  }

print.zeroyields <- function(x, ...)
  {
    cat("This is a dataset of zero-coupon yields.\n")
    cat(paste("Maturities range from", min(x$maturities), "to", max(x$maturities),"years.\n"))
    cat(paste("There are",nrow(x$yields), "observations between",x$dates[1], "and",x$dates[length(x$dates)],".\n"))
  }

summary.zeroyields <- function(object, ...)
  {
    print(summary(object$yields))
  }

plot.zeroyields <- function(x,...)
  {
    Z <- as.matrix(x$yields)
    X <- 1:nrow(x$yields)
    Y <- x$maturities

    open3d()
    persp3d(X, Y, Z, col = "green3", box = FALSE,xlab = "Dates", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
  }

print.dyntermstrc_yields <- function(x, ...){
  cat("---------------------------------------------------\n")
  cat("Parameters for yield based dynamic term structure estimation:\n")
  cat("---------------------------------------------------\n")
  cat("Method:",switch(x$method,"dl"="Diebold/Li","ns"="Nelson/Siegel","sv"="Svensson"),"\n")
  cat("Number of oberservations:",length(x$opt_result),"\n")
  cat("---------------------------------------------------\n")
  cat("Parameter summary:\n")
  cat("---------------------------------------------------\n")
  tsparam <- param.dyntermstrc_yields(x)
  print(lapply(tsparam,summary.default))
  cat("---------------------------------------------------\n")
}



summary.dyntermstrc_yields <- function(object, ...){
  x <- object
  y_mrsme <-  mean(sqrt(apply((x$yields-x$yhat)^2,2,mean)))
  y_maabse <- mean(apply(abs(x$yields-x$yhat),2,mean))
  sumry <- list()
  sumry$gof <- rbind(y_mrsme,y_maabse)
  rownames(sumry$gof) <- c("mean RMSE-Yields", "mean AABSE-Yields")

  ## extract convergence info
  for (i in (1:length(x$opt_result))) {
    sumry$convergence[i] <- x$opt_result[[i]]$convergence
    # TODO: solver message
  }
    
  class(sumry) <- "summary.dyntermstrc_yields"
  sumry
}


print.summary.dyntermstrc_yields <- function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")

    print.default(x$gof)

    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Convergence information from optim ():\n")
    cat("---------------------------------------------------\n")
    
    print.default(x$convergence)

  }

plot.dyntermstrc_yields <- function(x,...)
  {
    ## plot estimated yield curves in 3D
    sptrtfct <- switch(x$method,
                       "ns" = spr_ns,
                       "sv" = spr_sv)
    Z <- matrix(nrow=nrow(x$optparam),ncol=length(x$maturities))# OK
    for (i in 1:nrow(x$optparam)){
      Z[i,] <- sptrtfct(x$optparam[i,],x$maturities)
    }

    X <- 1:nrow(Z)
    Y <- x$maturities
    
    open3d()
    persp3d(X, Y, Z, col = "green3", box = FALSE,xlab = "Dates", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
    
  }


