
param <- function(obj,...) UseMethod("param") 


param.dyntermstrc <- function(x) {
  param <- t(mapply(function(i) x[[i]]$opt_result[[1]]$par, seq(length(x))))

  colnames(param) <- switch(x[[1]]$method,
                            "ns" = c("beta0","beta1","beta2","tau1"),
                            "sv" = c("beta0","beta1","beta2","tau1","beta3","tau2"),
                            "dl" = c("beta0","beta1","beta2"))
                           
                                       
  class(param) <- "dyntermstrc_param"
  param
}



print.dyntermstrc_param <- function(x,...){

  summary.default(x)

}


summary.dyntermstrc_param <- function(object,type="none",lags=1,selectlags="Fixed", ...) {
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
  
  class(sumry) <- "summary.dyntermstrc_param"
  sumry
}


print.summary.dyntermstrc_param <- function(x, ...) {
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


plot.dyntermstrc_param <- function(x,type="param",...){
  old.par <- par(no.readonly = TRUE) 
  param <- x
  
  # 2D plot of parameters
  if(type=="param") {
    if(ncol(x)==3) mfrow = c(1,3)
    if(ncol(x)==4) mfrow = c(2,2)
    if(ncol(x)==6) mfrow = c(2,3)

    par(mfrow=mfrow,...)
   
    plot(param[,1],type="l",xlab="Time",ylab=expression(hat(beta)[0]),
                col=1,lwd=2,... )
           grid()
           plot(param[,2],type="l",xlab="Time",ylab=expression(hat(beta)[1]),
           col=2,lwd=2,... )
           grid()
           plot(param[,3],type="l",xlab="Time",ylab=expression(hat(beta)[2]),
           col=3,lwd=2,... )
           grid()
    
    if(ncol(x)==4) {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
    }
    
    if(ncol(x)==6) {
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

 # 2D plot of parameter differences
  if(type=="diffparam") {
    if(ncol(x)==3) mfrow = c(1,3)
    if(ncol(x)==4) mfrow = c(2,2)
    if(ncol(x)==6) mfrow = c(2,3)

    par(mfrow=mfrow,...)

    diffparam <- apply(param,2,diff)

    for(i in seq(ncol(diffparam))) {
      plot(diffparam[,i],type="l",xlab="Time",
           ylab=colnames(diffparam)[i],col=i,lwd=2,... )
      grid()
    }
  }

    # ACF/PCF
  if(type=="acf") {
    if(ncol(x)==3) mfrow = c(2,3)
    if(ncol(x)==4) mfrow = c(4,2)
    if(ncol(x)==6) mfrow = c(4,3)

    par(mfrow=mfrow,...)
    
    if(ncol(x) > 3 ){
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
