
###################################################################
#               Dynamic term structure estimation                 #
###################################################################



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




###################################################################
#              print-method for summary.dyntermstrc               #
###################################################################






###################################################################
#                plot-method for dyntermstrc                      #
###################################################################



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









