###################################################################
#                        Spot rates cubic slines                  #
###################################################################

rm(list = ls())

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
T <- c(0, (T[i+1]-T[i-1])*((2*T[i+1-T[i]-T[i-1])/6+(t-T[i+1])/2)
       apply(as.matrix(m[,h]),2,max)
       +theta*(apply(as.matrix(m[,h+1]),2,max)-apply(as.matrix(m[,h]),2,max)),
       max(m[,ncol(m)]))


gi <- function(t,T,i,s){
  if(i==1){
    if(T[i]<=t&t<T[i+1]){
    g <- (T[i]-T[i-1])^2/6 + ((T[i]-T[i-1])*(t-T[i]))/2 + (t-T[i])^2/2 - (t-T[i])^3/(6*(T[i+1]-T[i]))}
  }
}
  
   # if(t>=T[i+1]){
   # print("Hello")
                                        #g <- (T[i+1]-T[i-1])#*((2*T[i+1]-T[i]-T[i-1])/6 + (t-T[i+1])/2)
    # }  
