rm(list = ls())

library(termstrc)
data(govbonds)
bdata <- prepro_bond("GERMANY",govbonds,c(0,10))
beta <- rep(1,4)
tau <- c(1,10)
m <- bdata$m[[1]]
cf <- bdata$cf[[1]]
w <- bdata$duration[[1]][,3]
p <- bdata$p[[1]]




objfct_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
      bsv <- c(beta[1:3],tau[1],beta[4],tau[2])
      phat <- bond_prices("sv",bsv,m,cf)$bond_prices
      loss_function(p, phat, w)
    }

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

objfct_sv_bonds_grid(beta, tau, m, cf, w, p)

control = list()
outer.iterations = 30
outer.eps = 1e-04

    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)

 lsparam <- constrOptim(theta = rep(0.01,4),
                                     f = objfct_sv_bonds_grid,
                                     #grad = grad_sv_bonds_grid,
                                     grad = NULL,
                                     ui = ui,
                                     ci = ci,
                                     mu = 1e-04,
                                     control = control,
                                     method = "Nelder-Mead",
                                     outer.iterations = outer.iterations,
                                     outer.eps = outer.eps,
                                     tau, m, cf, w, p) ## additional inputs for f and grad



 lsparam2 <- constrOptim(theta = rep(0.01,4),
                                     f = objfct_sv_bonds_grid,
                                     grad = grad_sv_bonds_grid,
                                     ui = ui,
                                     ci = ci,
                                     mu = 1e-04,
                                     control = control,
                                     method = "BFGS",
                                     outer.iterations = outer.iterations,
                                     outer.eps = outer.eps,
                                     tau, m, cf, w, p) ## additional inputs for f and grad


## ab hier schas

test1 <- function() {
  for (i in 1:10000) {
    objfct_sv_bonds_grid(beta, tau, m, cf, w, p)
  }
}

objfct_sv_bonds_grid2 <- function(beta, tau, m, cf, w, p) {
      bsv <- c(beta[1:3],tau[1],beta[4],tau[2])
      spot_rates <- spr_sv(bsv, m)/100
      spot_rates[is.nan(spot_rates)] <- 0        
      
      phat <- apply(cf*exp(-m*spot_rates), 2, sum)
      sum(w*((p - phat)^2))
    }

test2 <- function() {
  for (i in 1:10000) {
    objfct_sv_bonds_grid2(beta, tau, m, cf, w, p)
  }
}

system.time(test1())
system.time(test2())









Rprof()
test1()
Rprof(NULL)
summaryRprof()

Rprof()
test2()
Rprof(NULL)
summaryRprof()


