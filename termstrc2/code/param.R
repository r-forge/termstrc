
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
