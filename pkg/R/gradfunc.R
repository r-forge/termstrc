
# beta <- rep(1,4)
# tau <- rep(1,2)
# m <- matrix(1:10,10,10)
# cf <- matrix(c(rep(3,9),103),10,10)
# w <- rep(1,10)
# p <- rep(100,10)




grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p){

  a <- exp((-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m)*m/100)

  b <- -2*w*(p-apply(a*cf,2,sum))
  d <- a*cf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-apply(dm,2,sum)))
  gbeta2 <- sum(b*(-apply(d*tau[1]*(1-exp(-m/tau[1])),2,sum)))
  gbeta3 <- sum(b*(-apply(dm*(-exp(-m/tau[1]) +tau[1]*(1-exp(-m/tau[1]))/m),2,sum)))
  gbeta5 <- sum(b*(-apply(dm*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m),2,sum)))            


  c(gbeta1,gbeta2,gbeta3,gbeta5)
}


grad_sv_bonds <- function(beta,tau,m,cf,w,p){

   a <- exp((-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m)*m/100)

  b <- -2*w*(p-apply(a*cf,2,sum))
  d <- a*cf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-apply(dm,2,sum)))
  gbeta2 <- sum(b*(-apply(d*tau[1]*(1-exp(-m/tau[1])),2,sum)))
  gbeta3 <- sum(b*(-apply(dm*(-exp(-m/tau[1]) +tau[1]*(1-exp(-m/tau[1]))/m),2,sum)))     
  gbeta4 <- sum(b*(-apply(dm*(-beta[2]*exp(-m/tau[1])/tau[1] + beta[2]*(1-exp(-m/tau[1]))/m + beta[3]*(-exp(-m/tau[1])/tau[1]+ (1-exp(-m/tau[1]))/m - exp(-m/tau[1])*m/tau[1]^2)),2,sum)))
  gbeta5 <- sum(b*(-apply(dm*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m),2,sum)))   
  gbeta6 <- sum(b*(-apply(dm*beta[4]*( -exp(-m/tau[2])/tau[2] + (1-exp(-m/tau[2]))/m - exp(-m/tau[2])*m/tau[2]^2),2,sum)))

  c(gbeta1,gbeta2,gbeta3,gbeta4,gbeta5,gbeta6)
  
}

