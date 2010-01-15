summary.dyntermstrc <- function(object, ...) {
  x <- object

  # extract convergence info
  sumry <- list()
  sumry$convergence <- t(mapply(function(i) summary(x[[i]])$convergencegroup, seq_along(x)))
  colnames(sumry$convergence) <- x[[1]]$group 
  sumry$solvermsg <- t(mapply(function(i) summary(x[[i]])$convergence, seq_along(x)))
  colnames(sumry$solvermsg) <- x[[1]]$group 
  sumry$convprobs <- apply(sumry$convergence,2,function(x) which(x != "converged"))

  # adapt for multiple countries
  perrors <- list()
  yerrors <- list()
  
  perrors <- mapply(function(j) t(mapply(function(i) x[[i]]$perrors[[j]][,2], seq(length(x))))
  yerrors <- t(mapply(function(i) x[[i]]$yerrors[[j]][,2], seq(length(x))))
  p_mrsme <- mean(sqrt(apply(perrors^2,2,mean)))
  p_maabse <- mean(sqrt(apply(abs(perrors),2,mean)))
  y_mrsme <- mean(sqrt(apply(yerrors^2,2,mean)))
  y_maabse <- mean(sqrt(apply(abs(yerrors),2,mean)))
  
  sumry$gof <- rbind(p_mrsme,p_maabse,y_mrsme,y_maabse)
  colnames(sumry$gof) <- x[[1]]$group
  rownames(sumry$gof) <- c("mean RMSE-Prices", "mean AABSE-Prices", "mean RMSE-Yields", "mean AABSE-Yields")
    
  class(sumry) <- "summary.dyntermstrc"
  sumry
}



print.summary.dyntermstrc <- function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")

    print.default(x$gof)

    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Convergence information:\n")
    cat("---------------------------------------------------\n")
    
    if(length(x$convprobs)==0) {
      print.default("No convergence problems!")
    } else {
      cat("Convergence problems are at:\n")
      print.default(x$convprobs)
    }
   

 
    
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

print.dyntermstrc <- function(x,...){
  cat("---------------------------------------------------\n")
  cat("Parameters for dynamic term structure estimation:\n")
  cat("---------------------------------------------------\n")
  cat("Group: ",x[[1]]$group,"\n")
  cat("Method:",switch(x[[1]]$method,"dl"="Diebold/Li","ns"="Nelson/Siegel","sv"="Svensson"),"\n")
  cat("Number of oberservations:",length(x),"\n")
  cat("---------------------------------------------------\n")
  cat("Parameter summary:\n")
  cat("---------------------------------------------------\n")
  tsparam <- param.dyntermstrc(x)
  print(lapply(tsparam,summary.default))
  cat("---------------------------------------------------\n")
}
