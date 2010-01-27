## to 'lm' or not to 'lm'

y <- rnorm(100)
X <- cbind(rnorm(100),rnorm(100))
N <- 1000

system.time(for (i in 1:N) lm(y~X - 1)$coef)           
system.time(for (i in 1:N) solve(t(X)%*%X)%*%t(X)%*%y)  
