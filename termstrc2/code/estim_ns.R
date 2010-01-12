
estim_ns <- function(bonddata,group, matrange="all", method="ns", fit = "prices", weights="duration", startparam="auto",
           otype="nlminb",
           control=list(eval.max=1000,iter.max= 500),
           lambda=0.0609
           
           ) {

  # check inputs  
  if(fit=="yields"&weights!="none"){
  warning("For minimization of yield errors no weights are needed")
  weights <- "none"}

  # data preprocessing 

  prepro <- prepro_bond(group=group,bonddata=bonddata,
           matrange=matrange)


  n_group=prepro$n_group
  sgroup=prepro$sgroup
  positions=prepro$positions
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  yc=prepro$yc
  duration=prepro$duration
  durationc=prepro$durationc
 
  # automatically select start param
  if(!is.matrix(startparam)) {
   
  startparam <- t(mapply(function(i) c(mean(y[[i]][,2]),y[[i]][1,2] - mean(y[[1]][,2]),0,1,0,1) ,sgroup))
  
  rownames(startparam) <- group
  colnames(startparam) <- c("beta0","beta1","beta2","tau1","beta3","tau2")
  startparam <- switch(method,
                       "ns" = startparam[,1:4],
                       "sv"=startparam[,1:6],
                       "dl"=startparam[,1:3])
                         
  startparam <- matrix(startparam,nrow=length(sgroup))
    } else startparam <- startparam
  
  
  if(otype=="nlminb"){
  # objective function 
  obj_fct_prices <- function(b) {    # price error minimization
     loss_function(p[[k]],
     	bond_prices(method,b,m[[k]],cf[[k]],lambda)$bond_prices,duration[[k]][,3],weights)}
  
  obj_fct_yields <- function(b) {  
    loss_function(y[[k]][,2],bond_yields(rbind(
    -bond_prices(method,b,m[[k]],cf[[k]],bond_prices)$bond_prices,cf[[k]]),m_p[[k]])[,2],duration[[k]][,3],weights)} 
  
  
  obj_fct <- switch(fit,
                "prices" = obj_fct_prices,
                "yields" = obj_fct_yields)
                
  
  # lower and upper bounds for estimation   
  lower_bounds <- switch(method,
                       "ns" = c(0, -Inf, -Inf, 0),
                       "sv" = c(0, -Inf, -Inf, 0, -Inf, 0),
                       "dl" = c(0,-Inf,-Inf))
 
  upper_bounds <- switch(method,
                       "ns" = rep(Inf, 4),
                       "sv" = rep(Inf, 6),
                       "dl" = rep(Inf,3))
  } 
  # calculate optimal parameter vector
  opt_result <- list()

  if(otype=="nlminb") {
  for (k in sgroup){
    opt_result[[k]] = nlminb(startparam[k,],obj_fct,lower=lower_bounds, upper=upper_bounds,control=control)
                  
  }
}
     

 # data post processing 

  postpro <- postpro_bond(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,method,lambda)
 
  t=postpro$t
  
 
 # return list of results 
 result <- list(group=group,           # e.g. countries, rating classes
                 matrange=matrange,    # maturity range of bonds
                 method=method,        # method (Nelson/Siegel or Svensson)
                 fit=fit,              # fitting method (prices or yields)
                 weights=weights,      # weighting type for estimation
                 otype=otype,          # used optimisation algorithm
                 startparam=startparam,#calculated startparameters
                 n_group=n_group,      # number of groups,
                 spot=postpro$zcy_curves,      # zero coupon yield curves
                 spread=postpro$s_curves,      # spread curves
                 forward=postpro$fwr_curves,   # forward rate curves
                 discount=postpro$df_curves,   # discount factor curves
                 expoints=postpro$expoints,    # extrapolation points
       		 cf=cf,                # cashflow matrix
                 m=m,                  # maturity matrix
                 duration=duration,    # duration, modified duration, weights
                 durationc=durationc,  # duration based on clean yield
                 p=p,                  # dirty prices        
                 phat=postpro$phat,    # estimated dirty prices         
                 perrors=postpro$perrors,      # price errors
                 ac=ac,                # accrued interest
                 y=y,                  # maturities and yields
                 yc=yc,                # yield based on clean prices
                 yhat=postpro$yhat,            # estimated yields
                 yerrors=postpro$yerrors,      # yield errors
                 opt_result=opt_result, # optimisation results
                 lambda=lambda               
                 )
              
  # assign names to results list 
  for ( i in 9:length(result)) names(result[[i]]) <- group
    
  class(result) <- "termstrc_ns"
  result
 }
