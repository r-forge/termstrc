
###################################################################
#               Dynamic term structure estimation                 #
###################################################################

dyntermstrc <- function(dynbonddata,matrange,method,fit,weights,
                        startparam="auto",lambda=0.0609,otype="nlminb",bpeq="dirty",...) {

res <- list()
 
 
  # perform sequence of term structure estimations
  for (i in seq(length(dynbonddata))) {
    if(i>1){
      # use optimal parameters from previous period as start parameters
      b <- switch(method,
                  "Nelson/Siegel" = matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=4,byrow=TRUE),
                  "Svensson" =matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=6,byrow=TRUE),
                  "Diebold" =matrix(res[[i-1]]$opt_result[[1]]$par,nrow=1,ncol=3,byrow=TRUE))
      rownames(b) <- group
      colnames(b) <- switch(method,
                            "Nelson/Siegel" = c("beta0","beta1","beta2","tau1"),
                            "Svensson" = c("beta0","beta1","beta2","tau1","beta3","tau2"),
                            "Diebold" = c("beta0","beta2","beta3"))
                            
    } else b <- "auto"

    # static estimation
    group <- names(dynbonddata)[i]
    bonddata <- list()
    bonddata[[group]] <- dynbonddata[[i]]
    res[[i]] <- nelson_estim(group, bonddata=bonddata, matrange, 
                           method=method, fit, weights, startparam=b,
                           lambda=lambda,otype=otype,bpeq=bpeq,...)
  }
  class(res) <- "dyntermstrc"

  res
}

###################################################################
#               Parameters extractor function                     #
###################################################################

# TO DO:include time stamp

parameters <- function(x) {
  param <- t(mapply(function(i) x[[i]]$opt_result[[1]]$par, seq(length(x))))

  colnames(param) <- switch(x[[1]]$method,
                            "Nelson/Siegel" = c("beta0","beta1","beta2","tau1"),
                            "Svensson" = c("beta0","beta1","beta2","tau1","beta3","tau2"),
                            "Diebold" = c("beta0","beta1","beta2"))
                           
                                       
  class(param) <- "param"
  param
}

###################################################################
#          summary-method for dyntermstrc parameters              #
###################################################################

summary.param <- function(object,type="none",lags=1,selectlags="Fixed", ...) {
  x <- object

  sumry <- list()
  sumry$adflevels <- list()
  # Augmented Dickey Fuller Test for levels
  sumry$adflevels <- apply(x,2,function(x) ur.df(x,type=type,lags=lags,selectlags=selectlags)) #alternatively use adf.testx  
  sumry$adflevelsm <- matrix(NA,nrow=switch(type,"none"=3,"trend"=7),ncol=ncol(x))

  for (i in 1:length(sumry$adflevels)) {
    sumry$adflevelsm[switch(type,"none"=1,"trend"=c(1:3)),i] <- sumry$adflevels[[i]]@teststat # adf.test : $statisic
    sumry$adflevelsm[switch(type,"none"=2,"trend"=4),i] <- sumry$adflevels[[i]]@lags # adf.test: $parameter
    sumry$adflevelsm[switch(type,"none"=3,"trend"=c(5:7)),i] <- sumry$adflevels[[i]]@cval[,3] # adf.test: $p.value
  }
  rownames(sumry$adflevelsm) <- switch(type,"trend"=c(rep("Test statistic",3), "Lag order", "p-value-5pct","p-value-5pct","p-value-5pct"),"none"=c("Test statistic", "Lag order", "p-value-5pct"))
  colnames(sumry$adflevelsm) <- colnames(x)
  
  
  # Augmented Dickey Fuller Test for first differences
  sumry$adfdiff <- apply(x,2,function(x) ur.df(diff(x),type=type,lags=lags,selectlags=selectlags))
  sumry$adfdiffm <- matrix(NA,nrow=switch(type,"none"=3,"trend"=7),ncol=ncol(x))
  for (i in 1:length(sumry$adflevels)) {
    sumry$adfdiffm[switch(type,"none"=1,"trend"=c(1:3)),i] <- sumry$adfdiff[[i]]@teststat # adf.test : $statisic
    sumry$adfdiffm[switch(type,"none"=2,"trend"=4),i] <- sumry$adfdiff[[i]]@lags # adf.test: $parameter
    sumry$adfdiffm[switch(type,"none"=3,"trend"=c(5:7)),i] <- sumry$adfdiff[[i]]@cval[,3] # adf.test: $p.value
  }
  rownames(sumry$adfdiffm) <- switch(type,"trend"=c(rep("Test statistic",3), "Lag order", "p-value-5pct","p-value-5pct","p-value-5pct"),"none"=c("Test statistic", "Lag order", "p-value-5pct"))
  colnames(sumry$adfdiffm) <- colnames(x)

  sumry$paramcor <- cor(x)
  sumry$diffparamcor <- cor(apply(x,2,diff))
  
  class(sumry) <- "summary.param"
  sumry
}

###################################################################
#              print-method for summary.param                     #
###################################################################

print.summary.param <- function(x, ...) {
  cat("---------------------------------------------------\n")
  cat("ADF:\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  # Augmented Dickey Fuller Test for levels
  print.default(t(x$adflevelsm))
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("ADF of differences:\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  # Augmented Dickey Fuller Test for first differences
  print.default(t(x$adfdiffm))
  cat("\n")
  # correlation matrix of parameters
  cat("---------------------------------------------------\n")
  cat("Correlation of parameters:\n")
  cat("---------------------------------------------------\n")
  print.default(x$paramcor)
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("Correlation of differences:\n")
  cat("---------------------------------------------------\n")
  print.default(x$diffparamcor)
  cat("\n")

}


###################################################################
#              summary-method for dyntermstrc                     #
###################################################################

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


###################################################################
#              print-method for summary.dyntermstrc               #
###################################################################

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
}





###################################################################
#                plot-method for dyntermstrc                      #
###################################################################

plot.dyntermstrc <- function(x,type="param",mfrow = c(2,3),range=c(0,20), ...) {

  old.par <- par(no.readonly = TRUE) 
  
  # 2D plot of parameters
  if(type=="param") {
    if(x[[1]]$method=="Diebold") mfrow = c(1,3)
    if(x[[1]]$method=="Nelson/Siegel") mfrow = c(2,2)
    if(x[[1]]$method=="Svensson") mfrow = c(2,3)

    par(mfrow=mfrow,...)
    param <- parameters(x)

    #for(i in seq(ncol(param))) {
    #  plot(param[,i],type="l",xlab="Time",ylab=colnames(param)[i],
    #       col=i,lwd=2,... )
    #  grid()
    #}

    plot(param[,1],type="l",xlab="Time",ylab=expression(hat(beta)[0]),
                col=1,lwd=2,... )
           grid()
           plot(param[,2],type="l",xlab="Time",ylab=expression(hat(beta)[1]),
           col=2,lwd=2,... )
           grid()
           plot(param[,3],type="l",xlab="Time",ylab=expression(hat(beta)[2]),
           col=3,lwd=2,... )
           grid()
    
    if(x[[1]]$method=="Nelson/Siegel") {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
    }
    
    if(x[[1]]$method=="Svensson") {
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
    param <- parameters(x)

    X <- seq(if(range[1]==0) range[1]+0.1 else range[1],range[2],0.1)
    Y <- seq(nrow(param))

    Z <-  mapply(function(i) spotrates(method=x[[1]]$method,param[i,],X), seq(nrow(param)))

    persp(X,Y,Z,theta = -35, phi = 30, expand = 0.6, col = "lightgreen",
          ltheta = 120, shade = 0.55, ticktype = "detailed",xlab="Maturity",
          zlab="Zero-coupon yields",ylab="Time",box=TRUE,border=NA)

  }

  # 2D plot of parameter differences
  if(type=="diffparam") {
    if(x[[1]]$method=="Diebold") mfrow = c(1,3)
    if(x[[1]]$method=="Nelson/Siegel") mfrow = c(2,2)
    if(x[[1]]$method=="Svensson") mfrow = c(2,3)

    par(mfrow=mfrow,...)

    diffparam <- apply(parameters(x),2,diff)

    for(i in seq(ncol(diffparam))) {
      plot(diffparam[,i],type="l",xlab="Time",
           ylab=colnames(diffparam)[i],col=i,lwd=2,... )
      grid()
    }
  }

    # ACF/PCF
  if(type=="acf") {
    if(x[[1]]$method=="Diebold") mfrow = c(2,3)
    if(x[[1]]$method=="Nelson/Siegel") mfrow = c(4,2)
    if(x[[1]]$method=="Svensson") mfrow = c(4,3)

    par(mfrow=mfrow,...)
    
    #diffparam <- apply(parameters(x),2,diff)
    param <- parameters(x)
    if(x[[1]]$method %in% c("Svensson","Nelson/Siegel")){
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

###################################################################
#              Pricing errors extractor function                  #
###################################################################

perrors <- function(x) {
  # pricing errors
  errors <- t(mapply(function(i) x[[i]]$perrors[[1]][,2], seq(length(x))))
  # maturities at start date (needed for spacing plot method)
  ttmat <- x[[length(x)]]$y[[1]][,1]         
  perrors <- list(errors=errors,ttmat=ttmat)   
  class(perrors) <- "3derror"
  perrors
}


###################################################################
#                plot-method for 3derror                          #
###################################################################

# TODO: switch for yield errors

plot.3derror <- function(x, ... ) {
  old.par <- par(no.readonly = TRUE) 

  # errors <- switch(type,
  #                  "price" = perrors(x),
  #                  "yield" = yerrors(x))

  if (is.list(x)) {             # estimation error as list
    persp(seq(length(x$ttmat)), # use "x$ttmat" for maturity spacing => problems if bonds with same maturity
          seq(nrow(x$errors)),
          t(x$errors),
          theta = -35,
          phi = 30, expand = 0.6, col = "lightgreen",
          ltheta = 120, shade = 0.55, ticktype = "detailed",
          xlab="Maturity",zlab= " Price error",ylab="Time",
          box=TRUE,border=NA)
  } else                        # forecast error as matrix        
    persp(seq(ncol(x)),seq(nrow(x)),t(x), theta = -35, phi = 30,
          expand = 0.6, col = "lightgreen", ltheta = 120,
          shade = 0.55, ticktype = "detailed", xlab="Bond",
          zlab= " Price error",ylab="Time",box=TRUE,border=NA)

  on.exit(par(old.par))  
}


 
###################################################################
#              Yields extractor function                          #
###################################################################

yields <- function(x) {
  yields <- t(mapply(function(i) x[[i]]$y[[1]][,2], seq(length(x))))                  
  yields
}









