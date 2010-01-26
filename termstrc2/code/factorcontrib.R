factor_contrib_nss <- function(tparam,m){

fc_nss <- function(beta, m){
   fc1 <- beta[1]
   fc2 <- beta[2]*((1-exp(-m/beta[4]))/(m/beta[4]))
   fc3 <- beta[3]*(((1-exp(-m/beta[4]))/(m/beta[4]))-exp(-m/beta[4]))
   if(length(beta)==6) fc4 = beta[5]*(((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6])) else fc4=NULL
   return(list(fc1=fc1,fc2=fc2,fc3=fc3,fc4=fc4))
}

fc2 <- t(mapply(function(i) fc_nss(tparam[[1]][i,],m)$fc2, seq(nrow(tparam[[1]]))))
fc3 <- t(mapply(function(i) fc_nss(tparam[[1]][i,],m)$fc3, seq(nrow(tparam[[1]]))))


    Y <- seq(ncol(fc2))
    X <- seq(nrow(fc2))
    Z <- fc2
    open3d()
      persp3d(X,Y,Z,col = "green3",xlab="Time", zlab="Factor 2 -  Conbribution",ylab="Maturity (years)",box=FALSE)


     Y <- seq(ncol(fc3))
    X <- seq(nrow(fc3))
    Z <- fc3
    open3d()
      persp3d(X,Y,Z,col = "green3",xlab="Time", zlab="Factor3-  Conbribution",ylab="Maturity (years)",box=FALSE)

  if(ncol(tparam[[1]])==6) { fc4 <- t(mapply(function(i) fc_nss(tparam[[1]][i,],m)$fc4, seq(nrow(tparam[[1]])))) 
     Y <- seq(ncol(fc3))
     X <- seq(nrow(fc3))
     Z <- fc4
     open3d()
     persp3d(X,Y,Z,col = "green3",xlab="Time", zlab="Factor4-  Conbribution",ylab="Maturity (years)",box=FALSE)

                                                     
   Y <- seq(ncol(fc3))
    X <- seq(nrow(fc3))
    Z <- fc2 + fc4
    open3d()
      persp3d(X,Y,Z,col = "green3",xlab="Time", zlab="Factor 2+4 -  Conbribution",ylab="Maturity (years)",box=FALSE)
                           }
}










  
