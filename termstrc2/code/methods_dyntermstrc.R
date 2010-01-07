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




plot.dyntermstrc <- function(x,type="param",mfrow = c(2,3),range=c(0,20), ...) {

  old.par <- par(no.readonly = TRUE) 
  
  # 2D plot of parameters
  if(type=="param") {
    if(x[[1]]$method=="dl") mfrow = c(1,3)
    if(x[[1]]$method=="ns") mfrow = c(2,2)
    if(x[[1]]$method=="sv") mfrow = c(2,3)

    par(mfrow=mfrow,...)
    param <- param.dyntermstrc(x)

    plot(param[,1],type="l",xlab="Time",ylab=expression(hat(beta)[0]),
                col=1,lwd=2,... )
           grid()
           plot(param[,2],type="l",xlab="Time",ylab=expression(hat(beta)[1]),
           col=2,lwd=2,... )
           grid()
           plot(param[,3],type="l",xlab="Time",ylab=expression(hat(beta)[2]),
           col=3,lwd=2,... )
           grid()
    
    if(x[[1]]$method=="ns") {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
    }
    
    if(x[[1]]$method=="sv") {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
           plot(param[,5],type="l",xlab="Time",ylab=expression(hat(beta)[3]),
           col=5,lwd=2,... )
           grid()
           plot(param[,6],type="l",xlab="Time",ylab=expression(hat(tau)[2]),
           col=6,lwd=2,... )
           grid()
    }
    
  }

  # 3D plot of zero-coupon yield curves
  if(type=="3D") {
    param <- param.dyntermstrc(x)

    X <- seq(if(range[1]==0) range[1]+0.1 else range[1],range[2],0.1)
    Y <- seq(nrow(param))

    Z <-  mapply(function(i) spotrates(method=x[[1]]$method,param[i,],X,x[[1]]$lambda), seq(nrow(param)))

    persp(X,Y,Z,theta = -35, phi = 30, expand = 0.6, col = "lightgreen",
          ltheta = 120, shade = 0.55, ticktype = "detailed",xlab="Maturity",
          zlab="Zero-coupon yields",ylab="Time",box=TRUE,border=NA)

  }

  # 2D plot of parameter differences
  if(type=="diffparam") {
    if(x[[1]]$method=="dl") mfrow = c(1,3)
    if(x[[1]]$method=="ns") mfrow = c(2,2)
    if(x[[1]]$method=="sv") mfrow = c(2,3)

    par(mfrow=mfrow,...)

    diffparam <- apply(param.dyntermstrc(x),2,diff)

    for(i in seq(ncol(diffparam))) {
      plot(diffparam[,i],type="l",xlab="Time",
           ylab=colnames(diffparam)[i],col=i,lwd=2,... )
      grid()
    }
  }

    # ACF/PCF
  if(type=="acf") {
    if(x[[1]]$method=="dl") mfrow = c(2,3)
    if(x[[1]]$method=="ns") mfrow = c(4,2)
    if(x[[1]]$method=="sv") mfrow = c(4,3)

    par(mfrow=mfrow,...)
    
    param <- param.dyntermstrc(x)
    if(x[[1]]$method %in% c("sv","ns")){
     for(i in 1:(ncol(param)/2)) acf(param[,i],main=colnames(param)[i])
     for(i in 1:(ncol(param)/2)) pacf(param[,i],main=colnames(param)[i])
    
     for(i in (ncol(param)/2+ 1):ncol(param)) acf(param[,i],main=colnames(param)[i])
     for(i in (ncol(param)/2+ 1):ncol(param)) pacf(param[,i],main=colnames(param)[i])
    } else {

     for(i in 1:ncol(param)) acf(param[,i],main=colnames(param)[i])
     for(i in 1:ncol(param)) pacf(param[,i],main=colnames(param)[i])

    }
  }

on.exit(par(old.par))
}
