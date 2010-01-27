## Spot rates functions

## Nelson/Siegel spot rate function
spr_ns <- function(beta, m){
  (beta[1] + beta[2]*((1-exp(-m/beta[4]))/(m/beta[4])) +
   beta[3]*(((1-exp(-m/beta[4]))/(m/beta[4]))-exp(-m/beta[4])))
}

## Svensson spot rate function 
spr_sv <- function(beta, m){
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
   beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
   beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6])))
}

## Adjusted Svensson spot rate function
spr_asv <- function(beta, m){
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
   beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
   beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-(2*m)/beta[6])))
}

## Diebold/Li spot rate function
spr_dl <- function(beta,m,lambda){
  (beta[1] + beta[2]*((1-exp(-m*lambda))/(m*lambda))+
   beta[3]*(((1-exp(-m*lambda))/(m*lambda))-exp(-m*lambda)))
}

## Spot rate wrapper function
spotrates <- function(method,beta,m,lambda = 0.0609*12){ 
  switch(method,
         "ns" = spr_ns(beta,m),
         "sv" = spr_sv(beta,m),
         "asv"= spr_asv(beta,m),
         "dl" = spr_dl(beta,m,lambda))
}


## Forward rate functions

## Nelson/Siegel forward rate function
fwr_ns <- function(beta,m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
   beta[3]*(m/beta[4]*exp(-m/beta[4])))
}

## Svensson forward rate function
fwr_sv <- function(beta, m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
   beta[3] *m/beta[4]*exp(-m/beta[4]) +
   beta[5] *m/beta[6]*exp(-m/beta[6]))
}

## Adjusted Svensson forward rate function
fwr_asv <- function(beta, m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
   beta[3] *m/beta[4]*exp(-m/beta[4]) +
   beta[5] *(exp(-m/beta[6])+(2*m/beta[6] -1)*exp(-2*m/beta[6])))
}

## Diebold/Li forward rate function
fwr_dl <- function(beta, m,lambda) {
  (beta[1] + beta[2]*exp(-m*lambda)
   + beta[3]*(m*lambda*exp(-m*lambda)))
}

## Forward rates wrapper function
forwardrates <- function(method,beta,m,lambda){
  switch(method,
         "ns" = fwr_ns(beta,m),
         "sv" = fwr_sv(beta,m),
         "asv"= fwr_asv(beta,m),
         "dl"= fwr_dl(beta,m,lambda))
}


## Implied foreward rates calculation
impl_fwr <- function(m,s) {
  impl_fwr <- c(s[1],(s[-1]*m[-1] - s[-length(s)]*m[-length(m)])/(diff(m)))
  impl_fwr[1] <- impl_fwr[2]
  impl_fwr	
}

get_paramnames <- function(method){
  names <- c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")
  switch(method,"ns"= names[1:4],"sv"=names,"asv"=names,"dl"=names[1:3])
}

get_realnames <- function(method){
  switch(method,"dl"="Diebold/Li","ns"="Nelson/Siegel","sv"="Svensson","asv"="Adjusted Svensson")
}

### Loss function for parametric methods
get_objfct <- function(method) {
  objfct <- switch(method,
                   "ns" = objfct_ns,
                   "sv" = objfct_sv,
                   "asv" = objfct_asv)
  # TODO: replace with swich for spot rate function and remove objfct_*
}

### Gradient of loss function for parametric methods
get_grad_objfct <- function(method) {
  grad_objfct <- switch(method,
                   "ns" = grad_objfct_ns,
                   "sv" = grad_objfct_sv,
                   "asv" = grad_objfct_asv)
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
        sqrt(sum((y - spr_sv(beta,m))^2))
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

### Adjusted Svensson loss function for yields
objfct_asv <- function(beta, m, y)
      {
        sum((y - spr_asv(beta,m))^2)
      }

### Gradient of adjusted Svensson loss function for yields
grad_objfct_asv <- function(beta, m, y)
      {
        c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
          beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
          (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
          beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
          (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),
          
          sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
          
          sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
                 beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),

          sum(-2*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m)*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
          
          sum(-2*beta[5]*(-(1/(beta[6]*exp(m/beta[6]))) + (1 - exp(-m/beta[6]))/m - 
                          (2*m)/(beta[6]^2*exp((2*m)/beta[6])))*
              (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
               beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
               (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
          )
      }

### Constraints for constrOptim()

get_constraints <- function(method) {

  ## Diebold/Li
  
  if (method == "dl") {
    ui <- rbind(c(1,0,0),               # beta0 > 0
                c(1,1,0))               # beta0 + beta1 > 0
    ci <- c(0,0)
   }
  
  ## Nelson/Siegel
  
  if (method == "ns") {
    ui <- rbind(c(1,0,0,0),             # beta0 > 0
                c(1,1,0,0),             # beta0 + beta1 > 0
                c(0,0,0,1),             # tau1 > 0
                c(0,0,0,-1))            # tau1 < 30
    ci <- c(0,0,0,-30)
  }

  ## (Adjusted) Svensson

  if (method %in% c("sv","asv")) {
     ui <- rbind(c(1,0,0,0,0,0),        # beta0 > 0
                 c(1,1,0,0,0,0),        # beta0 + beta1 > 0
                 c(0,0,0,1,0,0),        # tau1 > 0
                 c(0,0,0,-1,0,0),       # tau1 < 30
                 c(0,0,0,0,0,1),        # tau2 > 0
                 c(0,0,0,0,0,-1),       # tau2 < 30
                 c(0,0,0,-1,0,1))       # tau2 - tau1 > 0
     ci <- c(0,0,0,-30,0,-30,0)
   }
    
  constraints <- list(ui = ui, ci = ci)
  constraints
}

