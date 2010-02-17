
# beta <- rep(1,4)
# tau <- rep(1,2)
# m <- matrix(1:10,10,10)
# cf <- matrix(c(rep(3,9),103),10,10)
# w <- rep(1,10)
# p <- rep(100,10)




grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p){

  a <- exp((-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m)*m/100)

  a[is.nan(a)] <- 0
  b <- -2*w*(p-colSums(a*cf))
  d <- a*cf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-colSums(dm)))
  gbeta2 <- sum(b*(-colSums(d*tau[1]*(1-exp(-m/tau[1])))))

  b3 <- dm*(-exp(-m/tau[1]) +tau[1]*(1-exp(-m/tau[1]))/m)
  b3[is.nan(b3)] <- 0
  gbeta3 <- sum(b*(-colSums(b3)))
  
  b5 <- dm*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)
  b5[is.nan(b5)] <- 0
  gbeta5 <- sum(b*(-colSums(b5)))            


  c(gbeta1,gbeta2,gbeta3,gbeta5)
}


grad_sv_bonds <- function(beta,tau,m,cf,w,p){

   a <- exp((-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m)*m/100)

  a[is.nan(a)] <- 0
  b <- -2*w*(p-apply(a*cf,2,sum))
  d <- a*cf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-apply(dm,2,sum)))
  gbeta2 <- sum(b*(-apply(d*tau[1]*(1-exp(-m/tau[1])),2,sum)))

  b3 <- dm*(-exp(-m/tau[1]) +tau[1]*(1-exp(-m/tau[1]))/m)
  b3[is.nan(b3)] <- 0
  gbeta3 <- sum(b*(-apply(b3,2,sum)))

  b4 <- dm*(-beta[2]*exp(-m/tau[1])/tau[1] + beta[2]*(1-exp(-m/tau[1]))/m + beta[3]*(-exp(-m/tau[1])/tau[1]+ (1-exp(-m/tau[1]))/m - exp(-m/tau[1])*m/tau[1]^2))
  b4[is.nan(b4)] <- 0
  gbeta4 <- sum(b*(-apply(b4,2,sum)))
   
  b5 <- dm*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)
  b5[is.nan(b5)] <- 0
  gbeta5 <- sum(b*(-apply(b5,2,sum)))   

  b6 <- dm*beta[4]*( -exp(-m/tau[2])/tau[2] + (1-exp(-m/tau[2]))/m - exp(-m/tau[2])*m/tau[2]^2)
  b6[is.nan(b6)] <- 0 
  gbeta6 <- sum(b*(-apply(b6,2,sum)))

  c(gbeta1,gbeta2,gbeta3,gbeta4,gbeta5,gbeta6)
  
}

