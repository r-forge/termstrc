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
  
  perrors <- mapply(function(j) t(mapply(function(i) x[[i]]$perrors[[j]][,2], seq(length(x)))),seq(x[[1]]$n_group),SIMPLIFY=FALSE)
  yerrors <- mapply(function(j) t(mapply(function(i) x[[i]]$yerrors[[j]][,2], seq(length(x)))),seq(x[[1]]$n_group),SIMPLIFY=FALSE)
  names(perrors) <- group
  names(yerrors) <- group
    
  p_mrsme <- mapply(function(i) mean(sqrt(apply(perrors[[i]]^2,2,mean))),seq_along(perrors))
  p_maabse <- mapply(function(i) mean(sqrt(apply(abs(perrors[[i]]),2,mean))),seq_along(perrors))
  y_mrsme <- mapply(function(i) mean(sqrt(apply(yerrors[[i]]^2,2,mean))),seq_along(yerrors))
  y_maabse <- mapply(function(i) mean(sqrt(apply(abs(yerrors[[i]]),2,mean))),seq_along(yerrors))
  
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
  
      
  # 3D plot of zero-coupon yield curves
    tsparam <- param.dyntermstrc(x)
    par(mfrow=c(1,nrow(tsparam[[1]])),...)
    X <- seq(if(range[1]==0) range[1]+0.1 else range[1],range[2],0.1)
    Y <- seq(nrow(tsparam[[1]]))

    for(j in seq_along(x)){
      Z <-  mapply(function(i) spotrates(method=x[[1]]$method,tsparam[[j]][i,],X,x[[1]]$lambda), seq(nrow(tsparam[[j]])))
      open3d()
      persp3d(X,Y,Z,col = "green3",xlab="Maturity (years)", zlab="Zero-coupon yields (in %)",ylab="Time",box=FALSE)
    }
 
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
