###################################################################
#                        Spot rates cubic slines                  #
###################################################################

rm(list = ls())
source("tools.R")

p <- c(96.6,93.71,91.56,90.24,89.74,90.04,91.09,92.82,95.19,98.14,101.60,105.54,109.90,114.64,119.73)

maturity <- seq(1,15)
coupon <- seq(2,9,0.5)

# cash flow matrix

N <- length(p)
cf <- matrix(0,N,N)
diag(cf) <- 100

for(i in 1:N) cf[1:maturity[i],i] <- cf[1:maturity[i],i]+coupon[i]

cf_p <- rbind(-p,cf)

# maturity matrix

m <- matrix(0,N,N)
diag(m) <- maturity

for(i in 1:N) m[1:i,i] <- seq(1,i)

m_p <- rbind(rep(0,N),m)

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


gi <- function(t,T,i,s){
  g <- rep(NA,length(t))
  for(j in 1:length(t)){
    if(i==1){
    if(T[i]<=t[j]&t[j]<T[i+1]){
     g[j] <- (T[i])^2/6 + ((T[i])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
    }
    if(t[j]>=T[i+1]){
     g[j] <- (T[i+1])*((2*T[i+1]-T[i])/6 + (t[j]-T[i+1])/2)
    }   
  }
  if(i>1&i<length(T)){
    if(t[j]<T[i-1]){
     g[j] <- 0
    }
    if(T[i-1]<=t[j]&t[j]<T[i]){
     g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
    }
    if(T[i]<=t[j]&t[j]<T[i+1]){
     g[j] <- (T[i]-T[i-1])^2/6 + ((T[i]-T[i-1])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
    }
    if(t[j]>=T[i+1]){
     g[j] <- (T[i+1]-T[i-1])*((2*T[i+1]-T[i]-T[i-1])/6 + (t[j]-T[i+1])/2)
    }
  }
   if(i==length(T)){
    if(t[j]<T[i-1]){
     g[j] <- 0
    }
    if(T[i-1]<=t[j]&t[j]<=T[i]){
     g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
    }
  } 
  if(i==s){
    g[j] <- t[j]
  }
}
  g
}


y <- apply(cf_p,2,sum)

X <- matrix(NA,N,s)

t = apply(m,2,max)

for(i in 1:s){
X[,i] <- apply(cf*gi(t,T,i,s),2,sum)
}

alpha <- coef(lm(-y~X-1))

t = seq(1,N,0.01)

dt <- rep(1,length(t))

for(i in 1:s){
  dt <- dt + alpha[i]*gi(t,T,i,s)
}

yhat <- -log(dt)/t


yields <- bond_yields(cf_p,m_p)
plot(yields[,1],yields[,2],ylim=c(0,0.08))
lines(t,yhat,type="l")
