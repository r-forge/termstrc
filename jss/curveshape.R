1m <- seq(0,70,0.01)     # time to maturity

beta_1 <- 1      # parameters
tau_1 <- 5
beta_2 <- 4
beta_0 <- 1
beta_3 <- 1
tau_2 <- 10



pdf("curveshape.pdf",width=9,height=9)
par(lwd=3)
# spot rates
plot(m,beta_0+ beta_1*((1-exp(-m/tau_1))/(m/tau_1))+ beta_2*(((1-exp(-m/tau_1))/(m/tau_1))-exp(-m/tau_1))+beta_3*(((1-exp(-m/tau_2))/(m/tau_2))-exp(-m/tau_2)),
     col=1,type="l",lty=1, ylim=c(0,3),
     xlab="Time to maturity",
     ylab="Model curves, spot rate")

# beta_0
lines(m,rep(beta_0,length(m)),col=2,lty=2)

# beta_1
lines(m, beta_1*((1-exp(-m/tau_1))/(m/tau_1)),type="l",col=3,lty=3)

# beta_2
lines(m, beta_2*(((1-exp(-m/tau_1))/(m/tau_1))-exp(-m/tau_1)),lty=4,col=4)

# beta_3
lines(m, beta_3*(((1-exp(-m/tau_2))/(m/tau_2))-exp(-m/tau_2)),lty=5,col=5)


legend("topright",legend=c(expression(s(m,b)),
       expression(beta[0]),
       expression(beta[1]*(frac(1-exp(-frac(m,tau[1])),frac(m,tau[1])))),
       expression(beta[2]*(frac(1-exp(-frac(m,tau[1])),frac(m,tau[1]))-exp(-frac(m,tau[1])))),
       expression(beta[3]*(frac(1-exp(-frac(m,tau[2])),frac(m,tau[2]))-exp(-frac(m,tau[2]))))
       ),
       lty=c(1,2,3,4,5), col=c(1,2,3,4,5),bty="n"
       )
dev.off()
