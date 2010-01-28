summary.dyntermstrc_nss <- function(object, ...) {
  x <- object

  # extract convergence info
  sumry <- list()
  sumry$convergence <- t(mapply(function(i) summary(x[[i]])$convergencegroup, seq_along(x)))
  rownames(sumry$convergence) <- x[[1]]$group 
  sumry$solvermsg <- t(mapply(function(i) summary(x[[i]])$convergence, seq_along(x)))
  rownames(sumry$solvermsg) <- x[[1]]$group 
 
  # adapt for multiple countries
  perrors <- list()
  yerrors <- list()
  
  perrors <- mapply(function(j) t(mapply(function(i) x[[i]]$perrors[[j]][,2], seq(length(x)))),seq(x[[1]]$n_group),SIMPLIFY=FALSE)
  yerrors <- mapply(function(j) t(mapply(function(i) x[[i]]$yerrors[[j]][,2], seq(length(x)))),seq(x[[1]]$n_group),SIMPLIFY=FALSE)
  names(perrors) <- x[[1]]$group
  names(yerrors) <- x[[1]]$group
    
  p_mrsme <- mapply(function(i) round(mean(sqrt(apply(perrors[[i]]^2,2,mean))),6),seq_along(perrors))
  p_maabse <- mapply(function(i) round(mean(apply(abs(perrors[[i]]),2,mean)),6),seq_along(perrors))
  y_mrsme <- mapply(function(i) round(mean(sqrt(apply((yerrors[[i]]*100)^2,2,mean))),6),seq_along(yerrors))
  y_maabse <- mapply(function(i) round(mean(apply(abs((yerrors[[i]]*100)^2),2,mean)),6),seq_along(yerrors))
  
  sumry$gof <- rbind(p_mrsme,p_maabse,y_mrsme,y_maabse)
  colnames(sumry$gof) <- x[[1]]$group
  rownames(sumry$gof) <- c("mean RMSE-Prices", "mean AABSE-Prices", "mean RMSE-Yields", "mean AABSE-Yields")
    
  class(sumry) <- "summary.dyntermstrc_nss"
  sumry
}



print.summary.dyntermstrc_nss <- function(x,...) {
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



plot.dyntermstrc_nss <- function(x,range=c(0,20), ...) {

  # 3D plot of zero-coupon yield curves
    tsparam <- param.dyntermstrc_nss(x)
    X <- seq(if(range[1]==0) range[1]+0.1 else range[1],range[2],0.1)
    Y <- seq(nrow(tsparam[[1]]))

    for(j in seq(x[[1]]$n_group)){
      Z <-  mapply(function(i) spotrates(method=x[[1]]$method,tsparam[[j]][i,],X,x[[1]]$lambda), seq(nrow(tsparam[[j]])))
      open3d()
      persp3d(X,Y,Z,col = "green3",xlab="Maturity (years)", zlab="Zero-coupon yields (in %)",ylab="Time",box=FALSE)
    }
 
}

print.dyntermstrc_nss <- function(x,...){
  cat("---------------------------------------------------\n")
  cat("Parameters for dynamic term structure estimation:\n")
  cat("---------------------------------------------------\n")
  cat("Group: ",x[[1]]$group,"\n")
  cat("Method:",get_realnames(x[[1]]$method),"\n")
  cat("Number of oberservations:",length(x),"\n")
  cat("---------------------------------------------------\n")
  cat("Parameter summary:\n")
  cat("---------------------------------------------------\n")
  tsparam <- param.dyntermstrc_nss(x)
  print(lapply(tsparam,summary.default))
  cat("---------------------------------------------------\n")
}
