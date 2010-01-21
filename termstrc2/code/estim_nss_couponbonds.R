##########################################################################
### Nelson/Siegel-type yield curve estimation method for 'couponbonds' ###
##########################################################################

estim_nss.couponbonds <- function(bonddata,                  # dataset (static)
                                  group,                     # names of countries for estimation c("Country 1", "Country 2", ...)
                                  matrange="all",            # maturity range in years c(min, max) 
                                  method="ns",
                                  startparam=NULL,           # startparameter matrix with columns c("beta0","beta1","beta2","tau1","beta3","tau2")
                                                             # otherwise globally optimal parameters are searched automatically
                                  lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                                  deltatau=0.1,              # interval for parameter grid
                                  #nlmbinbOptions = list(control = list()),
                                  constrOptimOptions = list(control = list(maxit = 2000), outer.iterations = 200, outer.eps = 1e-04)
           ) {

  ## data preprocessing
  prepro <- prepro_bond(group=group,bonddata=bonddata,matrange=matrange)

  n_group=prepro$n_group
  sgroup=prepro$sgroup
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  duration=prepro$duration
 
  ## automatically determine globally optimal start parameters
  spsearch <- list()
  length(spsearch) <- n_group
  
  if(is.null(startparam)){
    startparam <- matrix(ncol = 6, nrow = n_group)
    
    colnames(startparam) <- c("beta0","beta1","beta2","tau1","beta3","tau2")
    
    if (method == "dl") {startparam <- startparam[,1:3, drop=FALSE]}
    if (method == "ns") {startparam <- startparam[,1:4, drop=FALSE]}
    
    for (k in sgroup){
      print(paste("Searching startparameters for ", group[k]))
      spsearch[[k]] <- findstartparambonds(p[[k]],m[[k]],cf[[k]], duration[[k]][,3],
                                            method, deltatau)
      startparam[k,] <- spsearch[[k]]$startparam 
      print(startparam[k,])
    }
  }

  rownames(startparam) <- group
  
  ## objective function (weighted price error minimization) 
  obj_fct <- function(b) {
     loss_function(p[[k]],
     	bond_prices(method,b,m[[k]],cf[[k]],lambda)$bond_prices,duration[[k]][,3])}
                  
  ## calculate optimal parameter vectors
  opt_result <- list()

  for (k in sgroup){
    ## opt_result[[k]] <- constrOptim(theta = startparam[k,],
    ##                            f = obj_fct,
    ##                            grad = NULL,
    ##                            ui = ui,
    ##                            ci = ci,
    ##                            mu = 1e-04,
    ##                            control = control,
    ##                            method = "Nelder-Mead",
    ##                            outer.iterations = outer.iterations,
    ##                            outer.eps = outer.eps)
    opt_result[[k]] <- estimatezcyieldcurve(method, startparam[k,], obj_fct,constrOptimOptions) 
  }

  ## data post processing 
  postpro <- postpro_bond(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,method,lambda)
  
  ## return list of results 
  result <- list(group=group,                   # e.g. countries, rating classes
                 matrange=matrange,             # maturity range of bonds
                 method=method,                 # estimation method
                 startparam=startparam,         # calculated startparameters
                 n_group=n_group,               # number of groups,
                 lambda=lambda,                 # lambda parameter for dl
                 spsearch = spsearch,           # detailed data from start param search
                 spot=postpro$zcy_curves,       # zero coupon yield curves
                 spread=postpro$s_curves,       # spread curves
                 forward=postpro$fwr_curves,    # forward rate curves
                 discount=postpro$df_curves,    # discount factor curves
                 expoints=postpro$expoints,     # extrapolation points
       		 cf=cf,                         # cashflow matrix
                 m=m,                           # maturity matrix
                 duration=duration,             # duration, modified duration, weights
                 p=p,                           # dirty prices        
                 phat=postpro$phat,             # estimated dirty prices         
                 perrors=postpro$perrors,       # price errors
                 ac=ac,                         # accrued interest
                 y=y,                           # maturities and yields
                 yhat=postpro$yhat,             # estimated yields
                 yerrors=postpro$yerrors,       # yield errors
                 opt_result=opt_result          # optimisation results           
                 )
              
  for ( i in 6:length(result)) names(result[[i]]) <- group
  class(result) <- "termstrc_nss"
  result
}

### Estimate zero-coupon yield curve

estimatezcyieldcurve <- function(method, startparam, obj_fct, constrOptimOptions) {

  ## constraints

    if(method=="dl"){
    ui <- rbind(c(1,0,0),               # beta0 > 0
                c(1,1,0))               # beta0 + beta1 > 0
    ci <- c(0,0)
   }
 
  if(method=="ns"){
    ui <- rbind(c(1,0,0,0),             # beta0 > 0
                c(1,1,0,0),             # beta0 + beta1 > 0
                c(0,0,0,1),             # tau1 > 0
                c(0,0,0,-1))            # tau1 < 30
    ci <- c(0,0,0,-30)
    }

   if(method=="sv"){
     ui <- rbind(c(1,0,0,0,0,0),        # beta0 > 0
                 c(1,1,0,0,0,0),        # beta0 + beta1 > 0
                 c(0,0,0,1,0,0),        # tau1 > 0
                 c(0,0,0,-1,0,0),       # tau1 < 30
                 c(0,0,0,0,0,1),        # tau2 > 0
                 c(0,0,0,0,0,-1))       # tau2 < 30
     ci <- c(0,0,0,-30,0,-30)
    }


    ## lower <- switch(method,
    ##                 "ns" = c(0, -Inf, -Inf, 0),
    ##                 "sv" = c(0, -Inf, -Inf, 0, -Inf, 0),
    ##                 "dl" = c(0,-Inf,-Inf))
 
    ## upper <- switch(method,
    ##                 "ns" = rep(Inf, 4),
    ##                 "sv" = rep(Inf, 6),
    ##                 "dl" = rep(Inf,3))
  
    ## use nlminb() because performance is better
    ## opt_result <- nlminb(start = startparam,
    ##                    objective = obj_fct,
    ##                    grad = NULL,
    ##                    control = nlmbinbOptions$control,
    ##                    lower = lower,
    ##                    upper = upper)

    ##  use constrOptim() if b_0 + b_1 > 0 is not satisfied by nlminb()
    ##if(sum(opt_result$par[1:2])<0) {
   ##   warning("Constraint beta_0 + beta_1 > was violated by nlminb(), switching to constrOptim()")
      opt_result <- constrOptim(theta = startparam,
                                f = obj_fct,
                                grad = NULL,
                                ui = ui,
                                ci = ci,
                                mu = 1e-04,
                                control = constrOptimOptions$control,
                                method = "Nelder-Mead",
                                outer.iterations = constrOptimOptions$outer.iterations,
                                outer.eps = constrOptimOptions$outer.eps)
    #}
    opt_result
}

### Start parameter search routine for bond data

findstartparambonds <- function(p,m,cf, weights, method, deltatau = 0.1,
                                control = list(), outer.iterations = 30, outer.eps = 1e-04) {
  
  if(method=="dl"){
    startparam = rep(0.01,3)
    tau = NULL
    fmin = NULL
    optind = NULL
  }
 
  if(method=="ns"){
    tau <- seq(deltatau, max(m), deltatau) 
    fmin <- rep(NA, length(tau))
    lsbeta <- matrix(nrow = length(tau), ncol = 4)

    objfct <- function(b) {
      loss_function(p,bond_prices("dl",b,m,cf,1/tau[i])$bond_prices,weights)
    }

    ui <- rbind(c(1,0,0),                 # beta0 > 0
                c(1,1,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)
      
    for (i in 1:length(tau)){
      
      lsparam <- constrOptim(theta = rep(0.01,3), # start parameters for D/L, objective function is convex
                               f = objfct,
                               grad = NULL,
                               ui = ui,
                               ci = ci,
                               mu = 1e-04,
                               control = control,
                               method = "Nelder-Mead",
                               outer.iterations = outer.iterations,
                               outer.eps = outer.eps) 
      beta <- c(lsparam$par,tau[i])
      fmin[i] <- lsparam$value
      lsbeta[i,] <- beta 
    }
    optind <- which(fmin == min(fmin))
    startparam <- lsbeta[optind,]
  }
    
  if(method=="sv"){

    objfct <- function(b) {
      bsv <- c(b[1:3],tau1[i],b[4],tau2[j])
      loss_function(p,bond_prices("sv",bsv,m,cf)$bond_prices,weights)
    }

    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)
      
    tau1 <- seq(deltatau, max(m),deltatau)
    tau2 <- seq(deltatau, max(m),deltatau)
    tau <- cbind(tau1, tau2)
    fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
    lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
    for (i in 1:length(tau1))
      { 
        for (j in 1:length(tau2))
          {
            lsparam <- constrOptim(theta = rep(0.01,4),
                                   f = objfct,
                                   grad = NULL,
                                   ui = ui,
                                   ci = ci,
                                   mu = 1e-04,
                                   control = control,
                                   method = "Nelder-Mead",
                                   outer.iterations = outer.iterations,
                                   outer.eps = outer.eps)
              
            beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
            fmin[i,j] <- lsparam$value
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

### Startparameter grid search plots

plot.spsearch <- function(obj,...) {

  if(is.matrix(obj$tau)){
      contour(obj$tau[,1],obj$tau[,2],obj$fmin,nlevels=10,xlab = "tau_1", ylab = "tau_2",main = "Objective function")
      points(obj$tau[obj$optind[1],1],obj$tau[obj$optind[2],2],pch = 10, col = "red")
      open3d()
      persp3d(obj$tau[,1], obj$tau[,2], obj$fmin, col = "green3", box = FALSE,xlab = "tau_1", ylab = "tau_2", zlab = "Objective function")
      points3d(obj$tau[obj$optind[1],1],obj$tau[obj$optind[2],2],min(obj$fmin), col = "red")
  } else {
      plot(obj$tau,obj$fmin,xlab = "tau_1", ylab = "Objective function", type = "l")
      points(obj$tau[obj$optind],obj$fmin[obj$optind],pch = 10, col = "red")
  }
}
