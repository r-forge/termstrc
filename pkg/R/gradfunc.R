
grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p){

  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])
  t1emt1 <- tau[1]*(1 - emt1)
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-emt2 + tau[2]*(1 - emt2)/m)
  

  a <- exp((-beta[1] - beta[3]*emt1tm - beta[4]*emt2tm - (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))            


  c(gbeta1,gbeta2,gbeta3,gbeta5)
}


grad_sv_bonds <- function(beta,tau,m,cf,w,p){

  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])
  oemt1 <- (1-emt1)
  t1emt1 <- tau[1]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-emt2 + tau[2]*(1 - emt2)/m)
  emt1t <- emt1/tau[1]
  emt2t <- emt2/tau[2]

  a <- exp((-beta[1] - beta[3]*emt1tm - beta[4]*emt2tm - (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))     
  gbeta4 <- sum(b*(-cSums( dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/tau[1])),na.rm=TRUE)))
  gbeta6 <- sum(b*(-cSums(dm*beta[4]*( -emt2t + (1-emt2)/m - emt2t*m/tau[2]),na.rm=TRUE)))

  c(gbeta1,gbeta2,gbeta3,gbeta4,gbeta5,gbeta6)
  
}


cSums <- function (x, na.rm = FALSE, dims = 1L) {
    dn <- dim(x)
    n <- prod(dn[1L:dims])
    dn <- dn[-(1L:dims)]
    z <-  .Internal(colSums(x, n, prod(dn), na.rm))
    z
}

grad_ns_bonds_grid <- function(beta, tau, m, cf, w, p){
   emt1 <- exp(-m/tau[1])
   oemt1 <- (1-emt1)
   t1emt1 <- tau[1]*oemt1
   emt1tm <- (-emt1 + t1emt1/m)
   emt1t <- emt1/tau[1]
 
   a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
         (beta[2]*t1emt1)/m)*m)/100)

   acf <- a*cf
   b <- -2*w*(p-cSums(acf,na.rm=TRUE))
   d <- acf/100
   dm <- d*m

   gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
   gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
   gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))    

   c(gbeta1,gbeta2,gbeta3)
}


grad_ns_bonds <- function(beta, tau, m, cf, w, p){
   emt1 <- exp(-m/tau[1])
   oemt1 <- (1-emt1)
   t1emt1 <- tau[1]*oemt1
   emt1tm <- (-emt1 + t1emt1/m)
   emt1t <- emt1/tau[1]
 
   a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
         (beta[2]*t1emt1)/m)*m)/100)

   acf <- a*cf
   b <- -2*w*(p-cSums(acf,na.rm=TRUE))
   d <- acf/100
   dm <- d*m

   gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
   gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
   gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))    
   gbeta4 <- sum(b*(-cSums( dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/tau[1])),na.rm=TRUE)))
   
   c(gbeta1,gbeta2,gbeta3,gbeta4)
}



grad_asv_bonds_grid <- function(beta, tau, m, cf, w, p){
  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])
  oemt1 <- (1-emt1)
  t1emt1 <- tau[1]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-exp(-2*m/tau[2]) + (tau[2]*(1 -emt2))/m)
  emt1t <- emt1/tau[1]
 
  a <- exp(-(beta[1] + beta[3]*emt1tm + 
         beta[4]*emt2tm + 
         (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m

  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))     
  
  c(gbeta1, gbeta2, gbeta3, gbeta5)
  
}

grad_asv_bonds <-  function(beta, tau, m, cf, w, p){
  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])
  oemt1 <- (1-emt1)
  t1emt1 <- tau[1]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-exp(-2*m/tau[2]) + (tau[2]*(1 -emt2))/m)
  emt1t <- emt1/tau[1]
 
  a <- exp(-(beta[1] + beta[3]*emt1tm + 
         beta[4]*emt2tm + 
         (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m

  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1, na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta4 <- sum(b*(-cSums( dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/tau[1])),na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))     
  gbeta6 <- sum(b*(-cSums(dm*beta[4]*(-emt2/tau[2] + (1-emt2)/m - 2*exp(-2*m/tau[2])*m/tau[2]^2 ),na.rm=TRUE)))

  c(gbeta1, gbeta2, gbeta3, gbeta4, gbeta5, gbeta6)

}


