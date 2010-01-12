
estim_ns <- function(bonddata,group, matrange="all", method="ns", startparam="auto",
                     lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                     deltatau=0.1,              # interval for parameter grid
                     control=list(),            # options or optim() 
                     outer.iterations = 50,     # options for constrOptim()
                     outer.eps = 1e-05,
                     diagnosticplots = FALSE    # plots for start parameter search
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
  betastart <- matrix(ncol = 6, nrow = n_group)
  colnames(betastart) <- c("beta0","beta1","beta2","tau1","beta3","tau2")
  if (method == "dl") betastart <- betastart[,1:3]
  if (method == "ns") betastart <- betastart[,1:4]

  for (k in sgroup){
    print(paste("Searching startparameters for ", group[k]))
    betastart[k,] <- findstartparambonds(p[[k]],m[[k]],cf[[k]], duration[[k]][,3], method, deltatau, diagnosticplots, group[k])
    print(betastart[k,])
  }
  
  ## objective function (weighted price error minimization) 
  obj_fct <- function(b) {
     loss_function(p[[k]],
     	bond_prices(method,b,m[[k]],cf[[k]],lambda)$bond_prices,duration[[k]][,3])}
                  
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

  ## calculate optimal parameter vectors
  opt_result <- list()

  for (k in sgroup){
    opt_result[[k]] <- constrOptim(theta = betastart[k,],
                               f = obj_fct,
                               grad = NULL,
                               ui = ui,
                               ci = ci,
                               mu = 1e-04,
                               control = control,
                               method = "Nelder-Mead",
                               outer.iterations = outer.iterations,
                               outer.eps = outer.eps)              
  }

  ## data post processing 
  postpro <- postpro_bond(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,method,lambda)
  
  ## return list of results 
  result <- list(group=group,           # e.g. countries, rating classes
                 matrange=matrange,    # maturity range of bonds
                 method=method,        # method (Nelson/Siegel or Svensson)
                 startparam=startparam,#calculated startparameters
                 n_group=n_group,      # number of groups,
                 spot=postpro$zcy_curves,      # zero coupon yield curves
                 spread=postpro$s_curves,      # spread curves
                 forward=postpro$fwr_curves,   # forward rate curves
                 discount=postpro$df_curves,   # discount factor curves
                 expoints=postpro$expoints,    # extrapolation points
       		 cf=cf,                # cashflow matrix
                 m=m,                  # maturity matrix
                 duration=duration,    # duration, modified duration, weights
                 p=p,                  # dirty prices        
                 phat=postpro$phat,    # estimated dirty prices         
                 perrors=postpro$perrors,      # price errors
                 ac=ac,                # accrued interest
                 y=y,                  # maturities and yields
                 yhat=postpro$yhat,            # estimated yields
                 yerrors=postpro$yerrors,      # yield errors
                 opt_result=opt_result # optimisation results           
                 )
              
  for ( i in 6:length(result)) names(result[[i]]) <- group
  class(result) <- "termstrc_ns"
  result
}

findstartparambonds <- function(p,m,cf, weights, method, deltatau = 0.1,diagnosticplots = FALSE, name = "", control = list(), outer.iterations = 200, outer.eps = 1e-05) {
    if(method=="ns"){
      tau <- seq(deltatau,10,deltatau) # the first hump is within 10 years
      fmin <- rep(NA, length(tau))
      lsbeta <- matrix(nrow = length(tau), ncol = 4)

      ui <- rbind(c(1,0,0),             # beta0 > 0
            c(1,1,0))                   # beta0 + beta1 > 0
      ci <- c(0,0)

      objfct_dl <- function(b) {
        loss_function(p,
     	bond_prices("dl",b,m,cf,1/tau[i])$bond_prices,weights)
      }

      theta <- rep(0.01,3) # start parameters for D/L, should not matter because objective function is convex
      for (i in 1:length(tau)){
        lsparam <- constrOptim(theta = theta,
                               f = objfct_dl,
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

      if(diagnosticplots){
        X11()
        plot(tau,fmin,xlab = "tau_1", ylab = "Objective function", main = name, type = "l")
        points(tau[optind],fmin[optind],pch = 10, col = "red")
      }
    }
    
     if(method=="sv"){
     }
     startparam
}
