 
###################################################################
#            print-method for termstrc_singlecurve                #
###################################################################


print.termstrc_singlecurve <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for single-curve estimation:\n")
  cat("\n")
  cat("Method:",method,"\n")
  cat("Fitted:",fit,"\n")
  cat("Weights:",weights,"\n")
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  parameters <- mapply(function(i) x$opt_result[[i]]$par,1:length(x$opt_result))
  colnames(parameters) <- names(x$opt_result)
  n_par <- as.character(nrow(parameters))
  rownames(parameters) <- switch(n_par,
          "4"=c("beta_0","beta_1","beta_2","tau_1"),
          "6"=c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")) 
  print.default(parameters)
  cat("\n")
  }


 
###################################################################
#            plot-method for termstrc_singlecurve                 #
###################################################################

plot.termstrc_singlecurve <-
  function(x,maturity_min=1,maturity_max=15,spread_curves=TRUE,...) {
    # calculate theoretical yield curves
     yield_curves <- switch(x$method,
              "Nelson/Siegel" = mapply(function(i) 
              					nelson_siegel(x$opt_result[[i]]$par,
              					seq(maturity_min,maturity_max,0.01)),
              					1:x$n_countries),
              					
              "Svensson" = mapply(function(i) svensson(x$opt_result[[i]]$par,
              				seq(maturity_min,maturity_max,0.01)),1:x$n_countries))
              				
    # plot each yield curve seperately
    par(mfrow=c(2,2))
    for (i in 1:x$n_countries) {
      plot(seq(maturity_min,maturity_max,0.01),yield_curves[,i],
      type="l",
      ylim=c(0,0.05),
      xlab="Maturities",
      ylab="Yields",
      lwd=2,
      col="steelblue")
      title(names(x$opt_result)[i])
      grid()
      points(x$yields[[i]],col="red")
      par(ask=TRUE)
    }
    
    
    par(mfrow=c(1,2))
    
    ## plot all yield curves together
    matplot(seq(maturity_min,maturity_max,0.01),yield_curves,type="l",
    		col=1:x$n_countries,lty=1,lwd=2,
   
    xlab="Maturities",
    ylab="Yields",
    xlim=c(maturity_min,maturity_max))
    title("Yield curves")
    legend("bottomright",legend=names(x$opt_result),col=1:x$n_countries,lty=1,lwd=2)
    grid()
    
    ## calculate spread curves
    if (spread_curves==TRUE)
    yield_curve_ref <- yield_curves[,1]
    n <- x$n_countries-1              #number of spread curves
    spread_curves <- yield_curves[,2:(n+1)]-matrix(rep(yield_curve_ref,n),ncol=n)
    
    ### plot spread curves
    matplot(seq(maturity_min,maturity_max,0.01),spread_curves,type="l",col=2:(n+1),lty=1,lwd=2,
    xlab="Maturities",
    ylab="Spreads",
    xlim=c(maturity_min,maturity_max))
    title("Spread curves")
    legend("bottomright",legend=names(x$opt_result)[2:(n+1)],col=2:(n+1),lty=1,lwd=2)
    grid()
}  



###################################################################
#            summary-method for termstrc_singlecurve              #
###################################################################

summary.termstrc_singlecurve <-
    function(x) {
    RMSE_p <- mapply(function(i) rmse(x$dirty_prices[[i]],x$estimated_prices[[i]]),1:x$n_countries)
    AABSE_p <- mapply(function(i) aabse(x$dirty_prices[[i]],x$estimated_prices[[i]]),1:x$n_countries)
    RMSE_y <- mapply(function(i) rmse(x$yields[[i]][,2],x$spotrates[[i]]),1:x$n_countries)
    AABSE_y <- mapply(function(i) aabse(x$yields[[i]][,2],x$spotrates[[i]]),1:x$n_countries)
    
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$dirty_prices)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields","AABSE-Yields")
    cat("---------------------------------------------------\n")
    cat("Goodness of fit tests:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(gof)
    cat("\n")
    cat("\n")
    convergence <- as.matrix(mapply(function(i) x$opt_result[[i]]$message,1:length(x$opt_result)))
    colnames(convergence) <- "Convergence info"
    rownames(convergence) <- names(x$dirty_prices)
    print.default(convergence)
} 
