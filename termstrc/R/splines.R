###################################################################
#                        Cubic splines estimation                 #
###################################################################

splines_estim <-
  function(group,
           bonddata,
           maturity_spectrum="all"
           ) {

    
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
  cf_p <- mapply(function(k) create_cashflows_matrix(bonddata[[k]],include_price=TRUE),
                 1:n_group,SIMPLIFY=FALSE)
    
  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # create maturities matrix including zeros (needed for bond yield calculation)
  m_p <- mapply(function(k) create_maturities_matrix(bonddata[[k]],include_price=TRUE),
                1:n_group,SIMPLIFY=FALSE)
    
  # calculate dirty prices
  p <- mapply(function(k) bonddata[[k]]$PRICE + bonddata[[k]]$ACCRUED,1:n_group,SIMPLIFY=FALSE)
  
  # index for ordering
  positions <- mapply(function(k) order(apply(m[[k]],2,max)),1:n_group,SIMPLIFY=FALSE)
  
  # order matrices 
  cf <- mapply(function(k) cf[[k]][,positions[[k]]],1:n_group,SIMPLIFY=FALSE)
  cf_p <- mapply(function(k) cf_p[[k]][,positions[[k]]],1:n_group,SIMPLIFY=FALSE)
  m <- mapply(function(k) m[[k]][,positions[[k]]],1:n_group,SIMPLIFY=FALSE)
  m_p <- mapply(function(k) m_p[[k]][,positions[[k]]],1:n_group,SIMPLIFY=FALSE)
    
  # calculate bond yields	
  y <- mapply(function(k) bond_yields(cf_p[[k]],m_p[[k]]),
                   1:n_group,SIMPLIFY=FALSE)
    
  # Choosing knot points (McCulloch)
  K <- mapply(function(k) ncol(m[[k]]),1:n_group,SIMPLIFY=FALSE)
    
  # number of basis functions
  s <-  mapply(function(k) round(sqrt(K[[k]])),1:n_group,SIMPLIFY=FALSE)
    
  # only used for knot point finding
  i <- mapply(function(k) 2:(s[[k]]-2),1:n_group,SIMPLIFY=FALSE)  
  
  h <-  mapply(function(k) trunc(((i[[k]]-1)*K[[k]])/(s[[k]]-2)),1:n_group,SIMPLIFY=FALSE)
             
  theta <- mapply(function(k)((i[[k]]-1)*K[[k]])/(s[[k]]-2)-h[[k]],1:n_group,SIMPLIFY=FALSE)
    
  # knot points
  T <- mapply(function(k) c(0,
       apply(as.matrix(m[[k]][,h[[k]]]),2,max)
       +theta[[k]]*(apply(as.matrix(m[[k]][,h[[k]]+1]),2,max)-apply(as.matrix(m[[k]][,h[[k]]]),2,max)),
       max(m[[k]][,ncol(m[[k]])])),1:n_group,SIMPLIFY=FALSE)
 
  # use own knot points
  # T <- seq(0,15,1.5)
  # s <- length(T)+1
  
  # parameter estimation with OLS
  Y <- mapply(function(k) apply(cf_p[[k]],2,sum),1:n_group,SIMPLIFY=FALSE)
   
  X <- list()

  # k ... group index
  # j ... column index (bond)
  # sidx ... index for spline function  
 
  for (k in 1:n_group){  
  X[[k]] <- matrix(NA,ncol(m[[k]]),s[[k]])
                 
   for(sidx in 1:s[[k]]){
   X[[k]][,sidx] <- apply(cf[[k]]*mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]])),2,sum)
   }
  }
  
  alpha <- mapply(function(k) coef(lm(-Y[[k]]~X[[k]]-1)),1:n_group) # parameter vector
 
  # calculate discount factor matrix 
  dt <- list()
   for (k in 1:n_group){
   dt[[k]] <- matrix(1,nrow(m[[k]]),ncol(m[[k]]))
   for(sidx in 1:s[[k]]){
    dt[[k]] <- dt[[k]] + alpha[[k]][sidx]* mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]]))
   }
  }  

   browser()
  # calculate estimated prices 
  phat <- mapply(function(k) apply(cf[[k]]*dt[[k]],2,sum),1:n_group,SIMPLIFY=FALSE)
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-phat[[k]],cf[[k]]),m_p[[k]]),1:n_group,SIMPLIFY=FALSE)
  
 
  # calculate estimated zero coupon yield curves
  #t <- c(seq(0.01,0.4,0.1),seq(0.5,ceiling(max(mapply(function(i) max(y[[i]][,1]), 1:n_group))),0.01))
  
  t <- mapply(function(k) seq(0.01, max(T[[k]]),0.01), 1:n_group) 
  
  zcy_curves <- matrix(NA,nrow=length(t),ncol=n_group+1)
  zcy_curves[,1] <- t 
          
  dt_zcy <- list()
  for( k in 1:n_group) {
   dt_zcy[[k]] <- rep(1,length(t[[k]]))
    for(sidx in 1:s[[k]]){  
    dt_zcy[[k]] <- dt_zcy[[k]] + alpha[[k]][sidx]*gi(t,T[[k]],sidx,s[[l]])
   }
  zcy_curves[,l+1] <- -log(dt_zcy[[l]])/t          
  }
  
  # calculate spread curves              	    
 	if(n_group != 1) { 
   scurves <- zcy_curves[,3:(n_group+1)] - zcy_curves[,2] 	    
    } else scurves = "none" 
 
 #browser()
 # return list of results
 result <- list(  group=group,          # e.g. countries, rating classes
                  matrange=matrange,    # maturity range of bonds
                  n_group=n_group,      # number of groups,
                  T=T,                  #knot points 
                  zcy_curves=zcy_curves, # zero coupon yield curves
                  scurves=scurves,      # spread curves
                  cf=cf,                # cashflow matrix
                  m=m,                  # maturity matrix
                  p=p,                  # dirty prices
                  phat=phat,            # estimated prices
                  y=y,                  # maturities and yields
                  yhat=yhat,            # estimated yields
                  alpha=alpha           # cubic splines parameters                             
                 )
                 
  # assign names to results list 
  for ( i in 5:length(result)) names(result[[i]]) <- names(bonddata)
    
  class(result) <- "cubicsplines"
  result
 
}



