
## Calculation of the duration,modified duration and duration based weights 
duration <- function (cf_p,m_p,y){
  y <- matrix(rep(y,nrow(m_p)),ncol=ncol(m_p),byrow=TRUE)
  # mac cauly duration
  d <- apply(cf_p*m_p*exp(-y*m_p),2,sum)/-cf_p[1,]
  # modified duration
  md <- d/(1+y[1,])
  omega <- (1/d)/sum(1/d)
  dur <- cbind(d,md,omega)
  colnames(dur) <- c("Duration","Modified duration","Weights")
  dur
}

## Root mean squared error
rmse <- function (actual,estimated) {
  e <- actual - estimated
  sqrt(mean(e^2))
}

## Average absolute error
aabse <-function (actual,estimated){
  e <- actual - estimated	
  mean(abs(e))
}   							


## Loss function: mean squared (weighted error)
loss_function <- function(p,phat,omega) {
  sum(omega*((p-phat)^2))
}
