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

  a[is.nan(a)] <- 0
  b <- -2*w*(p-apply(a*cf,2,sum))
  d <- a*cf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-apply(dm,2,sum)))
  gbeta2 <- sum(b*(-apply(d*tau[1]*(1-exp(-m/tau[1])),2,sum)))

  b3 <- dm*(-exp(-m/tau[1]) +tau[1]*(1-exp(-m/tau[1]))/m)
  b3[is.nan(b3)] <- 0
  gbeta3 <- sum(b*(-apply(b3,2,sum)))
  
  b5 <- dm*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)
  b5[is.nan(b5)] <- 0
  gbeta5 <- sum(b*(-apply(b5,2,sum)))            


  c(gbeta1,gbeta2,gbeta3,gbeta5)
}

objfct_sv_bonds_grid(beta, tau, m, cf, w, p)

control = list()
outer.iterations = 30
outer.eps = 1e-04

ui <- rbind(c(1,0,0,0),                 # beta0 > 0
            c(1,1,0,0))                 # beta0 + beta1 > 0
ci <- c(0,0)

system.time(
lsparam <- constrOptim(theta = rep(0.01,4),
                       f = objfct_sv_bonds_grid,
                       grad = NULL,
                       ui = ui,
                       ci = ci,
                       mu = 1e-04,
                       control = control,
                       method = "Nelder-Mead",
                       outer.iterations = outer.iterations,
                       outer.eps = outer.eps,
                       tau, m, cf, w, p) ## additional inputs for f and grad
)

system.time(
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
)

## ab hier schas

test1 <- function() {
  for (i in 1:5) {
    lsparam <- constrOptim(theta = rep(0.01,4),
                       f = objfct_sv_bonds_grid,
                       grad = NULL,
                       ui = ui,
                       ci = ci,
                       mu = 1e-04,
                       control = control,
                       method = "Nelder-Mead",
                       outer.iterations = outer.iterations,
                       outer.eps = outer.eps,
                       tau, m, cf, w, p) ## additional inputs for f and grad
  }
}

test2 <- function() {
  for (i in 1:5) {
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
  }
}

system.time(test1())
system.time(test2())



objfct_sv_bonds_grid2 <- function(beta, tau, m, cf, w, p) {
      bsv <- c(beta[1:3],tau[1],beta[4],tau[2])
      spot_rates <- spr_sv(bsv, m)/100
      spot_rates[is.nan(spot_rates)] <- 0        
      
      phat <- apply(cf*exp(-m*spot_rates), 2, sum)
      sum(w*((p - phat)^2))
    }





Rprof()
test1()
Rprof(NULL)
summaryRprof()

Rprof()
test2()
Rprof(NULL)
summaryRprof()


