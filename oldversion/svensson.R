
###################################################################
#                        Spot Rates Nelson/Siegel                 #
###################################################################

nelson_siegel <-
  function(beta, m, curve="spotrate") {

  switch(curve,
  spotrate = beta[1] + beta[2] * ((1-exp(-m/beta[4]))/(m/beta[4]))
             + beta[3]*(((1-exp(-m/beta[4]))/(m/beta[4]))
             -exp(-m/beta[4])),
  forwardrate = beta[1]+beta[2]*exp(-m/beta[4])+beta[3]*(m/beta[4])
                *exp(-m/beta[4]),
  creditspread = beta[1] + beta[2] * ((1-exp(-m/beta[4]))/(m/beta[4]))
             + beta[3]*exp(-m/beta[4])             
 )                    
}

###################################################################
#                        Spot Rates Svensson                      #
###################################################################

svensson <-
  function(beta, m) {
  beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
  beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
  beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6]))}

###################################################################
#                        loss function                            #
###################################################################

loss_function <- function(p,phat,omega,weights) {
  if (weights=="none") omega <- rep(1,length(p))
  sum(omega*((p-phat)^2))}

###################################################################
#                 Term structure estimation                       #
###################################################################

# OPTIONS: - countries: vector with countries, i.e. c("AUSTRIA","GERMANY")
#						first country is the reference country for 
#						credit spread estimation
#          - bonddata: dataset where countries are included
#          - maturity_spectrum: include only bonds within that range
#                               "all" ... use all available bonds
#          - method: "Nelson/Siegel", "Svensson"
#          - fit: "prices", "yields"
#          - weights: "none", "duration"
#          - startparam: matrix with startparameters; one row for one
#			 country; number of columns is 4 for the nelson/siegel
#			 and 6 for the svensson method 		
#		   - control: a list with control parameters for the 
#			 nlminb optimisation function	   

termstrc_estim <-
  function(countries,
           bonddata,
           maturity_spectrum="all",
           method="Nelson/Siegel",
           fit = "prices",
           weights="none",
           startparam,
           control=list(eval.max=1000)
           
           ) {

  # select given countries from bonddata
  bonddata <- bonddata[countries]
   
  if (length(maturity_spectrum=1)){bonddata <- bonddata} 
  else{bonddata <- maturity_range(bonddata,maturity_spectrum[1],maturity_spectrum[2])}
  
  # number of countries 
  n_countries <- length(bonddata) 
    
  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)

  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # calculate dirty prices
  p <- mapply(function(i) bonddata[[i]]$PRICE + bonddata[[i]]$ACCRUED,1:n_countries)
  names(p) <- names(bonddata)   
 
  #?
  cf_p <- mapply(function(i) create_cashflows_matrix(bonddata[[i]],include_price=TRUE),
                 1:n_countries)
  names(cf_p) <- names(bonddata)
  #?
  m_p <- mapply(function(i) create_maturities_matrix(bonddata[[i]],include_price=TRUE),
                1:n_countries)
  names(m_p) <- names(bonddata)
  
  # calculate bond yields	
  yields <- mapply(function(i) bond_yields(cf_p[[i]],m_p[[i]]),
                   1:n_countries)
  names(yields) <- names(bonddata)
  
  #calculate duration
  duration <- mapply(function(i) duration(cf_p[[i]],m_p[[i]],yields[[i]][,2]),
                   1:n_countries)
  names(duration) <- names(bonddata)
  
  # SINGLE-CURVE ESTIMATION 
  # objective function type
 
  obj_fct_prices <- function(b) {
     loss_function(p[[i]],
     	bond_prices(method,b,m[[i]],cf[[i]])$bond_prices,duration[[i]][,3],weights)}
  
  obj_fct_yields <- function(b) {  
    loss_function(yields[[i]][,"Yield"],bond_yields(rbind(
    -bond_prices(method,b,m[[i]],cf[[i]])$bond_prices,cf[[i]]),m_p[[i]]))}
    
   obj_fct <- switch(fit,
                "prices" = obj_fct_prices,
                "yields" = obj_fct_yields)
  
 # bounds for single curve estimation   
 lower_bounds<-switch(method,
                       "Nelson/Siegel" = c(0, -Inf, -Inf, 0),
                       "Svensson" = c(0, -Inf, -Inf, 0, -Inf, 0))
 
 upper_bounds<-switch(method,
                       "Nelson/Siegel" = rep(Inf, 4),
                       "Svensson" = rep(Inf, 6))
 
 # single curve estimation
 #calculate optimal parameter 
 opt_result <- list()
 for (i in 1:n_countries){
 opt_result[[i]] <- nlminb(startparam[i,],obj_fct, 
 						lower = lower_bounds, upper = upper_bounds,control=control)
 }
 names(opt_result) <- names(bonddata)
 
 #estimated prices
 
 estimated_prices <- mapply(function(i) bond_prices(method,opt_result[[i]]$par,
                     			 m[[i]],cf[[i]])$bond_prices,1:n_countries)
 names(estimated_prices) <- names(bonddata)
 
 #calculate spotrates according to chosen approach
  
 spotrates <- mapply(function(i) srates(method,opt_result[[i]]$par,
 									yields[[i]][,1]),1:n_countries) 
 
 names(spotrates)<-names(bonddata)						
 #return list of results 
 result <- list(maturity_spectrum=maturity_spectrum,method=method,
       			fit=fit,weights=weights,n_countries=n_countries,
       			cashflows=cf,maturities=m,dirty_prices=p,
                estimated_prices=estimated_prices,yields=yields,
                opt_result=opt_result,spotrates=spotrates)
 
 class(result) <- "termstrc_singlecurve"    
 result


 }
 
 
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



###################################################################
#                   Bond pricing function                         #
###################################################################

bond_prices <-
  function(method="Nelson/Siegel", beta, m, cf,multicurve=FALSE) {
  
  
  
  # calculate spot rates for single-curve estimation
  
  spot_rates <-
    switch(method,
           "Nelson/Siegel" = nelson_siegel(beta,m),
           "Svensson"=svensson(beta,m)
           )
   # replace NAs by zeroes
   spot_rates[is.nan(spot_rates)] <- 0        
     
  
  # calculate discount factors
  discount_factors <- exp(-m*spot_rates)

  # calculate bond prices
  bond_prices <- apply(cf*discount_factors, 2, sum)  
  
  
  # return spot rates, discount factors and bond prices
  return (list(spot_rates=spot_rates,
               discount_factors=discount_factors,
               bond_prices=bond_prices))

}
    

	



