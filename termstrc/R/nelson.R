###################################################################
#               Nelson-Siegel, Svensson  estimation               #
###################################################################
   
nelson_estim <-
  function(group,
           bonddata,
           maturity_spectrum="all",
           method="Nelson/Siegel",
           fit = "prices",
           weights="none",
           startparam,
           control=list(eval.max=1000)
           
           ) {

  # check inputs  
  if(fit=="yields"&weights!="none"){
  warning("For minimization of yield errors no weights are needed")}
    
  # select given group from bonddata
  bonddata <- bonddata[group]
  
  # select data according to chosen maturity range
  if (length(maturity_spectrum)==1) {bonddata <- bonddata }else
   {bonddata <- maturity_range(bonddata,maturity_spectrum[1],maturity_spectrum[2]) }

  # number of groups 
  n_group <- length(bonddata) 
    
  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)

  # create cashflows matrix including dirty price (needed for bond yield calculation)
  cf_p <- mapply(function(i) create_cashflows_matrix(bonddata[[i]],include_price=TRUE),
                 1:n_group,SIMPLIFY=FALSE)

  names(cf_p) <- names(bonddata)

  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # create maturities matrix including zeros (needed for bond yield calculation)
  m_p <- mapply(function(i) create_maturities_matrix(bonddata[[i]],include_price=TRUE),
                1:n_group,SIMPLIFY=FALSE)
  names(m_p) <- names(bonddata)
  
  # calculate dirty prices
  p <- mapply(function(i) bonddata[[i]]$PRICE + bonddata[[i]]$ACCRUED,1:n_group,SIMPLIFY=FALSE)
  names(p) <- names(bonddata)   
 
  # calculate bond yields	
  yields <- mapply(function(i) bond_yields(cf_p[[i]],m_p[[i]]),
                   1:n_group,SIMPLIFY=FALSE)
  names(yields) <- names(bonddata)
  
  # calculate duration
      
  duration <- mapply(function(i) duration(cf_p[[i]],m_p[[i]],yields[[i]][,2]),
                   1:n_group,SIMPLIFY=FALSE)
  names(duration) <- names(bonddata)
  
  # objective function
 
  obj_fct_prices <- function(b) {    # price error minimization
     loss_function(p[[i]],
     	bond_prices(method,b,m[[i]],cf[[i]])$bond_prices,duration[[i]][,3],weights)}
  
  obj_fct_yields <- function(b) {    # yield error minimization
    loss_function(yields[[i]][,"Yield"],spotrates(method,b,yields[[i]][,"Maturity"]))}
    
  obj_fct <- switch(fit,
                "prices" = obj_fct_prices,
                "yields" = obj_fct_yields)
  
  # lower and upper bounds for estimation   
  lower_bounds <- switch(method,
                       "Nelson/Siegel" = c(0, -Inf, -Inf, 0),
                       "Svensson" = c(0, -Inf, -Inf, 0, -Inf, 0))
 
  upper_bounds <- switch(method,
                       "Nelson/Siegel" = rep(Inf, 4),
                       "Svensson" = rep(Inf, 6))
 
  # calculate optimal parameter vector
  opt_result <- list()

  for (i in 1:n_group){
    opt_result[[i]] <- nlminb(startparam[i,],obj_fct, 
     lower = lower_bounds, upper = upper_bounds,control=control)
  }
  names(opt_result) <- names(bonddata)
 
  # theoretical bond prices with estimated parameters
  estimated_prices <- mapply(function(i) bond_prices(method,opt_result[[i]]$par,
       m[[i]],cf[[i]])$bond_prices,1:n_group,SIMPLIFY=FALSE)
  names(estimated_prices) <- names(bonddata)
 
  # calculate spotrates according to chosen approach
  spot_rates <- mapply(function(i) spotrates(method,opt_result[[i]]$par,
 		    yields[[i]][,1]),1:n_group,SIMPLIFY=FALSE) 
 
  names(spot_rates)<-names(bonddata)						

  # return list of results 
  result <- list(maturity_spectrum=maturity_spectrum,method=method,
       		fit=fit,weights=weights,n_group=n_group,
       		cashflows=cf,maturities=m,dirty_prices=p,duration=duration,
                estimated_prices=estimated_prices,yields=yields,
                opt_result=opt_result,spot_rates=spot_rates)
 

  class(result) <- switch(method,
 	"Nelson/Siegel" = "nelson",
 	"Svensson" = "svensson")

  result
   }


###################################################################
#                        Spot rates Nelson/Siegel                 #
###################################################################

nelson_siegel <-
  function(beta, m) {
    (beta[1] + beta[2]*((1-exp(-m/beta[4]))/(m/beta[4]))
    + beta[3]*(((1-exp(-m/beta[4]))/(m/beta[4]))-exp(-m/beta[4])))}

###################################################################
#                        Spot rates Svensson                      #
###################################################################

svensson <-
  function(beta, m) {
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
  beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
  beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6])))}

###################################################################
#                        loss function                            #
###################################################################

loss_function <-
  function(p,phat,omega,weights="none") {
  if (weights=="none") omega <- rep(1,length(p))
  sum(omega*((p-phat)^2))}

	



