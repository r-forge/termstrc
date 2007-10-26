beta_1 <- 1             # parameters
tau_1 <- 1
beta_2 <- 1
beta_0 <- 1

#par(lwd=2)

m <- seq(0,10,0.01)

# beta_1
#plot(m,beta_1*exp(-m/tau_1),type="l",col=1,
#     xlab="Time to maturity",
#     ylab="Model curves",ylim=c(0,2.1),lty=4)

# beta_2
#lines(m,beta_2*(m/tau_1)*exp(-m/tau_1),lty=2,col=2)

# forward rate curve
#lines(m,beta_0+beta_1*exp(-m/tau_1)+beta_2*(m/tau_1)*exp(-m/tau_1),col=3)

#abline(h=1,lty=3,col=4)

#legend("topright",legend=c(expression(beta[1]*exp(-m/tau[1])),
#       expression(beta[2]*(m/tau[1])*exp(-m/tau[1])),
#       expression(beta[0]),
#       expression(f(m,b))),
#       lty=c(4,2,3,1), col=c(1,2,4,3))


nelson <- function(m,b){ b[1]+ b[2]*((1-exp(-m/b[4]))/(m/b[4]))+ b[3]*(((1-exp(-m/b[4]))/(m/b[4]))-exp(-m/b[4]))}

b <- cbind(5,-1,-12:12,1)


plot(m, nelson(m,b[1,]), ylim=c(0,8),type="l", xlab="Maturity", ylab="Spotrate curves")
for(i in 2:nrow(b)) lines(m,nelson(m,b[i,]))

svensson <-
  function(beta, m) {
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
  beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
  beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6])))}




#for(i in 2:nrow(b)) lines(m,svensson(b[i,],m))

par(mfrow=c(2,4))

b <- cbind(5,-1,1,1,4,4:0,1)

plot(m, svensson(b[1,],m),type="l",xlab="",ylab="")

plot(m,svensson(c(5,-1,1,1,1,1,10),m), type="l", xlab="",ylab="")

plot(m,svensson(c(5,-4,-15,1,0.1,6,4),m), type="l",xlab="",ylab="")

#plot(m,svensson(c(5,-4,-1,10,4,0.5,7),m), type="l")
plot(m,svensson(c(5,-4,-1,15,4,0.5,4),m), type="l",xlab="",ylab="")

plot(m,svensson(c(5,-4,-1,10,15,0.5,1),m), type="l",xlab="",ylab="")

plot(m,svensson(c(5,-4,-5,15,15,0.1,1),m), type="l",xlab="",ylab="")

plot(m,svensson(c(5,-4,-5,15,10,2,10),m), type="l",xlab="",ylab="")

plot(m,svensson(c(5,-4,1 ,1,4,2,1),m), type="l",xlab="",ylab="")

