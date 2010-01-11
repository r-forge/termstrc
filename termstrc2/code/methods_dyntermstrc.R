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

  perrors <- t(mapply(function(i) x[[i]]$perrors[[1]][,2], seq(length(x))))
  yerrors <- t(mapply(function(i) x[[i]]$yerrors[[1]][,2], seq(length(x))))
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
  cat("Method:",x[[1]]$method,"\n")
  cat("Fitted:",x[[1]]$fit,"\n")
  cat("Weights:",x[[1]]$weights,"\n")
  cat("Number of oberservations:",length(x),"\n")
  cat("Number of bonds:",ncol(x[[1]]$cf[[1]]),"\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  tsparam <- param.dyntermstrc(x)
  print(summary.default(tsparam))
  cat("---------------------------------------------------\n")
}
