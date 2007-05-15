
###################################################################
#                        Spot rates Nelson/Siegel                 #
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
#                        Spot rates Svensson                      #
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
#			 nlminb optimization function	   

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
  
  # select data according to chosen maturity range
  if (length(maturity_spectrum)==1) {bonddata <- bonddata }else
   {bonddata <- maturity_range(bonddata,maturity_spectrum[1],maturity_spectrum[2]) }


  # number of countries 
  n_countries <- length(bonddata) 
    
  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)

  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # calculate dirty prices
  p <- mapply(function(i) bonddata[[i]]$PRICE + bonddata[[i]]$ACCRUED,1:n_countries,SIMPLIFY=FALSE)
  names(p) <- names(bonddata)   
 
  #
  cf_p <- mapply(function(i) create_cashflows_matrix(bonddata[[i]],include_price=TRUE),
                 1:n_countries,SIMPLIFY=FALSE)
  names(cf_p) <- names(bonddata)
  #
  m_p <- mapply(function(i) create_maturities_matrix(bonddata[[i]],include_price=TRUE),
                1:n_countries,SIMPLIFY=FALSE)
  names(m_p) <- names(bonddata)
  
  # calculate bond yields	
  yields <- mapply(function(i) bond_yields(cf_p[[i]],m_p[[i]]),
                   1:n_countries,SIMPLIFY=FALSE)
  names(yields) <- names(bonddata)
  
  #calculate duration
  duration <- mapply(function(i) duration(cf_p[[i]],m_p[[i]],yields[[i]][,2]),
                   1:n_countries,SIMPLIFY=FALSE)
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
                     			 m[[i]],cf[[i]])$bond_prices,1:n_countries,SIMPLIFY=FALSE)
 names(estimated_prices) <- names(bonddata)
 
 #calculate spotrates according to chosen approach
  
 spotrates <- mapply(function(i) srates(method,opt_result[[i]]$par,
 		    yields[[i]][,1]),1:n_countries,SIMPLIFY=FALSE) 
 
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
    

	



