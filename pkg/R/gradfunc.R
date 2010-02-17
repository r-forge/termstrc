
# beta <- rep(1,4)
# tau <- rep(1,2)
# m <- matrix(1:10,10,10)
# cf <- matrix(c(rep(3,9),103),10,10)
# w <- rep(1,10)
# p <- rep(100,10)




grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p){

  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])

  a <- exp((-beta[1] - beta[3]*(-emt1 + (tau[1]*(1 - emt1))/m) - beta[4]*(-emt2 + (tau[2]*(1 - emt2))/m) - (beta[2]*tau[1]*(1 - emt1))/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-colSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-colSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-colSums(d*tau[1]*(1-emt1),na.rm=TRUE)))
  gbeta3 <- sum(b*(-colSums(dm*(-emt1 +tau[1]*(1-emt1)/m),na.rm=TRUE)))
  gbeta5 <- sum(b*(-colSums(dm*(-emt2 + (tau[2]*(1 - emt2))/m),na.rm=TRUE)))            


  c(gbeta1,gbeta2,gbeta3,gbeta5)
}


grad_sv_bonds <- function(beta,tau,m,cf,w,p){

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

  b4 <- dm*(-beta[2]*exp(-m/tau[1])/tau[1] + beta[2]*(1-exp(-m/tau[1]))/m + beta[3]*(-exp(-m/tau[1])/tau[1]+ (1-exp(-m/tau[1]))/m - exp(-m/tau[1])*m/tau[1]^2))
  b4[is.nan(b4)] <- 0
  gbeta4 <- sum(b*(-colSums(b4)))
   
  b5 <- dm*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)
  b5[is.nan(b5)] <- 0
  gbeta5 <- sum(b*(-colSums(b5)))   

  b6 <- dm*beta[4]*( -exp(-m/tau[2])/tau[2] + (1-exp(-m/tau[2]))/m - exp(-m/tau[2])*m/tau[2]^2)
  b6[is.nan(b6)] <- 0 
  gbeta6 <- sum(b*(-colSums(b6)))

  c(gbeta1,gbeta2,gbeta3,gbeta4,gbeta5,gbeta6)
  
}

