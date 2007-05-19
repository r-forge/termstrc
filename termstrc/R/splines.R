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

browser()
  
# Choosing knot points (McCulloch)

K <- N

# number of basis functions
s <- round(sqrt(K))

i = 2:(s-2)

h <- trunc(((i-1)*K)/(s-2))
theta <- ((i-1)*K)/(s-2)-h

# knot points
T <- c(0,
       apply(as.matrix(m[,h]),2,max)
       +theta*(apply(as.matrix(m[,h+1]),2,max)-apply(as.matrix(m[,h]),2,max)),
       max(m[,ncol(m)]))

# parameter estimation with OLS
y <- apply(cf_p,2,sum)
X <- matrix(NA,N,s)
t = apply(m,2,max)

for(i in 1:s){
X[,i] <- apply(cf*gi(t,T,i,s),2,sum)
}

alpha <- coef(lm(-y~X-1))   # parameter vector


dt <- rep(1,length(t))      # disccount factors

for(i in 1:s){
  dt <- dt + alpha[i]*gi(t,T,i,s)  
}

spot_rates <- -log(dt)/t           # estimated yields

estimated_prices <- apply(cf*dt,2,sum)
}


#yields <- bond_yields(cf_p,m_p)
#plot(yields[,1],yields[,2],ylim=c(0,0.08))
#lines(t,spot_rates,type="l")
