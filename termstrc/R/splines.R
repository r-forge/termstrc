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

  positions <- mapply(function(i) order(apply(m[[i]],2,max)),1:n_group,SIMPLIFY=FALSE)
  names(positions) <- names(bonddata)


  cf <- mapply(function(i) cf[[i]][,positions[[i]]],1:n_group,SIMPLIFY=FALSE)
  cf_p <- mapply(function(i) cf_p[[i]][,positions[[i]]],1:n_group,SIMPLIFY=FALSE)
  m <- mapply(function(i) m[[i]][,positions[[i]]],1:n_group,SIMPLIFY=FALSE)
  m_p <- mapply(function(i) m_p[[i]][,positions[[i]]],1:n_group,SIMPLIFY=FALSE)
 
  
  # calculate bond yields	
  yields <- mapply(function(i) bond_yields(cf_p[[i]],m_p[[i]]),
                   1:n_group,SIMPLIFY=FALSE)
  names(yields) <- names(bonddata)

  
# Choosing knot points (McCulloch)

K <- mapply(function(i) ncol(m[[i]]),1:n_group,SIMPLIFY=FALSE)
names(K) <- names(bonddata)
  
# number of basis functions
s <-  mapply(function(i) round(sqrt(K[[i]])),1:n_group,SIMPLIFY=FALSE)
names(s) <- names(bonddata)
  
i <- mapply(function(i) 2:(s[[i]]-2),1:n_group,SIMPLIFY=FALSE)  # only used for knot point finding
names(i) <- names(bonddata)

h <-  mapply(function(j) trunc(((i[[j]]-1)*K[[j]])/(s[[j]]-2)),1:n_group,SIMPLIFY=FALSE)
names(h) <- names(bonddata)
            
theta <- mapply(function(j)((i[[j]]-1)*K[[j]])/(s[[j]]-2)-h[[j]],1:n_group,SIMPLIFY=FALSE)
names(theta) <- names(bonddata)

  
# knot points
T <- mapply(function(i) c(0,
       apply(as.matrix(m[[i]][,h[[i]]]),2,max)
       +theta[[i]]*(apply(as.matrix(m[[i]][,h[[i]]+1]),2,max)-apply(as.matrix(m[[i]][,h[[i]]]),2,max)),
       max(m[[i]][,ncol(m[[i]])])),1:n_group,SIMPLIFY=FALSE)
names(T) <- names(bonddata)

# use own knot points
# T <- seq(0,15,1.5)
# s <- length(T)+1
  
# parameter estimation with OLS
y <- mapply(function(i) apply(cf_p[[i]],2,sum),1:n_group,SIMPLIFY=FALSE)
names(y) <- names(bonddata)

  
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

names(X) <- names(bonddata)
  
alpha <- mapply(function(k) coef(lm(-y[[k]]~X[[k]]-1)),1:n_group) # parameter vector

dt <- list()
  
for (k in 1:n_group){
  dt[[k]] <- matrix(1,nrow(m[[k]]),ncol(m[[k]]))
  for(sidx in 1:s[[k]]){
  dt[[k]] <- dt[[k]] + alpha[[k]][sidx]* mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]]))
  }
}  

browser()

k <-1
t <- seq(0,30,0.5)

dt2 <- list()
dt2[[k]] <- rep(1,length(t))
               
for(sidx in 1:s[[k]]){  
dt2[[k]] <- dt2[[k]] + alpha[[k]][sidx]*gi(t,T[[k]],sidx,s[[k]])
}

spot_rates <- -log(dt2[[k]])/t           # estimated yields

estimated_prices <- list()
estimated_prices <- mapply(function(k) apply(cf[[k]]*dt[[k]],2,sum),1:n_group,SIMPLIFY=FALSE)

yhat <- mapply(function(k) bond_yields(rbind(-estimated_prices[[k]],cf[[k]]),m_p[[k]]),1:n_group,SIMPLIFY=FALSE)


 k = 2 
 plot(yields[[k]],ylim=c(0,0.06))
 #lines(m[[k]][,ncol(m[[k]])],-log(dt[[k]][,ncol(m[[k]])])/m[[k]][,ncol(m[[k]])]) 
 points(yhat[[k]][,1],yhat[[k]][,2],col="blue")
 plot(t,spot_rates,type="l")
 points(yields[[k]][,1],yields[[k]][,2]) 
}

#-log(dt[[k]][,ncol(m[[k]])])/m[[k]][,ncol(m[[k]])]

#yields <- bond_yields(cf_p,m_p)
#plot(yields[,1],yields[,2],ylim=c(0,0.08))
#lines(t,spot_rates,type="l")
