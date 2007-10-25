beta_1 <- 1             # parameters
tau_1 <- 1
beta_2 <- 1
beta_0 <- 1

par(lwd=2)

m <- seq(0,10,0.01)

# beta_1
plot(m,beta_1*exp(-m/tau_1),type="l",col=1,
     xlab="Time to maturity",
     ylab="Model curves",ylim=c(0,2.1),lty=4)

# beta_2
lines(m,beta_2*(m/tau_1)*exp(-m/tau_1),lty=2,col=2)

# forward rate curve
lines(m,beta_0+beta_1*exp(-m/tau_1)+beta_2*(m/tau_1)*exp(-m/tau_1),col=3)

abline(h=1,lty=3,col=4)

legend("topright",legend=c(expression(beta[1]*exp(-m/tau[1])),
       expression(beta[2]*(m/tau[1])*exp(-m/tau[1])),
       expression(beta[0]),
       expression(f(m,b))),
       lty=c(4,2,3,1), col=c(1,2,4,3))


spotrate <- function(m,b){ b[1]+ b[2]*((1-exp(-m/b[4]))/(m/b[4]))+ b[3]*(((1-exp(-m/b[4]))/(m/b[4]))-exp(-m/b[4]))}

b <- c( 5,1,1,1)

b <- cbind(5,-1,-12:12,1)


plot(m,spotrate(m,b[1,]), ylim=c(-1,8),type="l")
for(i in 2:nrow(b)) lines(m,spotrate(m,b[i,]))
