### Nelson/Siegel and Svensson spot curve estimaton from zero yields 

zeroyields <- function(maturities, yields, dates)
  {
    zy <- list(maturities = maturities, yields = yields, dates = dates)
    class(zy) <- "zeroyields"
    zy
  }

print.zeroyields <- function(obj)
  {
    cat("This is a dataset of zero-coupon yields.\n")
    cat(paste("Maturities range from", min(obj$maturities), "to", max(obj$maturities),"years.\n"))
    cat(paste("There are",nrow(obj$yields), "observations between",obj$dates[1], "and",obj$dates[length(obj$dates)],".\n"))
  }

summary.zeroyields <- function(obj)
  {
    print(summary(obj$yields))
  }

plot.zeroyields <- function(obj)
  {
    z <- as.matrix(obj$yields)
    x <- 1:nrow(obj$yields)
    y <- obj$maturities

    open3d()
    persp3d(x, y, z, col = "green3", box = FALSE,xlab = "Dates", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
  }

estim_nss <- function(obj, ...) UseMethod("estim_nss")

estim_nss.zeroyields <- function (obj, method = "ns", deltatau = 1)
  {
    if(method=="ns"){
          spsearch <- findstartparamyields(obj$yields[1,],obj$maturities, method, deltatau)
          startparam <- spsearch$startparam
          objfct <- objfct_ns
          grad_objfct <- grad_objfct_ns
          ui <- rbind(c(1,0,0,0),  # beta0 > 0
                      c(1,1,0,0),  # beta0 + beta1 > 0
                      c(0,0,0,1),  # tau1 > 0
                      c(0,0,0,-1)) # tau1 < 30
          ci <- c(0,0,0,-30)
    }

    if(method=="sv"){
      
          spsearch <- findstartparamyields(obj$yields[1,],obj$maturities, method, deltatau)
          startparam <- spsearch$startparam
          
          objfct <- objfct_sv
          grad_objfct <- grad_objfct_sv
          ui <- rbind(c(1,0,0,0,0,0),  # beta0 > 0
                      c(1,1,0,0,0,0),  # beta0 + beta1 > 0
                      c(0,0,0,1,0,0),  # tau1 > 0
                      c(0,0,0,-1,0,0), # tau1 < 30
                      c(0,0,0,0,0,1),  # tau2 > 0
                      c(0,0,0,0,0,-1)) # tau2 < 30
          ci <- c(0,0,0,-30,0,-30)
    }

    optresult <- list()
    optparam <- matrix(nrow = nrow(obj$yields), ncol = length(startparam)) 
    for (i in 1:nrow(obj$yields)){
    
      yields <- obj$yields[i,]

      if (i<2)
        beta <- startparam
      else
        beta <- optparam[i-1,]
        
      optresult[[i]] <- estimateyieldcurve(yields, obj$maturities, beta, objfct, grad_objfct, ui, ci)
      optparam[i,] <- optresult[[i]]$par
    }
    colnames(optparam) <- switch(method,
          "ns"=c("beta_0","beta_1","beta_2","tau_1"),
          "sv"=c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")) 
    result <- list(optparam = optparam, optresult = optresult, method = method,
                   maturities = obj$maturities, dates = obj$dates, spsearch = spsearch, yields = obj$yields)
    class(result) <- "termstrc_yields"
    result
  }

plot.termstrc_yields <- function(obj)
  {
    ## plot parameters
    nparam <- ncol(obj$optparam) 
    op <- par(mfrow = c(2,nparam/2))
    for (i in 1:nparam){
      plot(obj$optparam[,i], type = "l", xlab = "Time", ylab = colnames(obj$optparam)[i],lwd=2,col=i)
      grid()
    }
    par(op)

    ## plot estimated yield curves in 3D
    sptrtfct <- switch(obj$method,
                       "ns" = spr_ns,
                       "sv" = spr_sv)
    z = matrix(nrow=nrow(obj$optparam),ncol=length(obj$maturities))
    for (i in 1:nrow(obj$optparam)){
      z[i,] <- sptrtfct(obj$optparam[i,],obj$maturities)
    }

    x <- 1:nrow(z)
    y <- obj$maturities
    
    open3d()
    persp3d(x, y, z, col = "green3", box = FALSE,xlab = "Dates", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
    
  }

estimateyieldcurve <- function(y, m, beta, objfct, grad_objfct, ui, ci)
  {    
    opt_result <- constrOptim(theta = beta,
                               f = objfct,
                               grad = grad_objfct,
                               ui = ui,
                               ci = ci,
                               mu = 1e-04,
                               control = list(),
                               method = "BFGS",
                               outer.iterations = 50,
                               outer.eps = 1e-05,
                               m,y)
          
  }

findstartparamyields <- function(y,m, method, deltatau = 0.1)
  {
    if(method=="ns"){
      tau <- seq(0.1,max(m),deltatau)
      fmin <- rep(NA, length(tau))
      lsbeta <- matrix(nrow = length(tau), ncol = 4)
      for (i in 1:length(tau)){
        X <- cbind(rep(1,length(y)),
                   ((1 - exp(-m/tau[i]))/(m/tau[i])),
                   (((1 - exp(-m/tau[i]))/(m/tau[i])) - exp(-m/tau[i])))

        lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
        beta <- c(lsparam[1:3],tau[i])
        fmin[i] <- objfct_ns(beta, m, y)
        lsbeta[i,] <- beta 
      }
      optind <- which(fmin == min(fmin))
      startparam <- lsbeta[optind,]
    }
    
     if(method=="sv"){
       tau1 <- seq(0.1, max(m), deltatau)
       tau2 <- seq(0.1, max(m), deltatau)
       tau <- cbind(tau1, tau2)
       fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
       lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
       for (i in 1:length(tau1))
         {
           for (j in 1:length(tau2))
             {
               ## reparametrize to avoid nonsingular matrix
               if (i == j){
                      X <- cbind(rep(1,length(y)),
                          ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                          - exp(-m/tau1[i]))

               lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
               beta <- c(lsparam[1],lsparam[2]-lsparam[3],lsparam[3]/2,tau1[i],lsparam[3]/2,tau2[j])
             } else
               {
                X <- cbind(rep(1,length(y)),
                          ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                          (((1 - exp(-m/tau1[i]))/(m/tau1[i])) - exp(-m/tau1[i])),
                          (((1 - exp(-m/tau2[j]))/(m/tau2[j])) - exp(-m/tau2[j])))

               lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
               beta <- c(lsparam[1:3],tau1[i],lsparam[4],tau2[j])
               }
               ## check parameter contraints (beta_0 > 0, beta_0 + beta_1 > 0, beta[2] should not explode)
               if(beta[1]>0 && ((beta[1]+beta[2])>0 && beta[2]<20)){
                 fmin[i,j] <- objfct_sv(beta, m, y)
             } else{
               fmin[i,j] <- 1
             }
               lsbeta[(i-1)*length(tau1)+j,] <- beta 
             } 
         }
       optind <- which(fmin == min(fmin),arr.ind=TRUE)
       startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]
     }
   
  result <- list(startparam = startparam, tau = tau, fmin = fmin, optind = optind)
  class(result) <- "spsearch"
  result
  }

### Nelson/Siegel loss function for yields
objfct_ns <- function(beta, m, y)
      {
        sum((y - spr_ns(beta,m))^2)
      }

### Gradient of Nelson/Siegel loss function for yields
grad_objfct_ns <- function(beta, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
          (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
           (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),

          sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
           (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
           (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
                 beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))
              *(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
                (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
          )
      }

### Svensson loss function for yields
objfct_sv <- function(beta, m, y)
      {
        sum((y - spr_sv(beta,m))^2)
      }

### Gradient of Svensson loss function for yields
grad_objfct_sv <- function(beta, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
        beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
        (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),

          sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
      beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(-2*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m)*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(-2*beta[5]*(-(1/(beta[6]*exp(m/beta[6]))) + (1 - exp(-m/beta[6]))/m - m/(beta[6]^2*exp(m/beta[6])))*
    (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
      beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
      (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
          )
      }











