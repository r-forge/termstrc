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
  
