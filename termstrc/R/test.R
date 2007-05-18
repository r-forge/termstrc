rm(list = ls())
load("eurobonds.RData")
source("nelson.R")
source("tools.R")

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
maturity_spectrum <- "all"
method <- "Nelson/Siegel"
fit <- "yields"
weights <- "none"
control <- list(eval.max=100000)

b <- matrix(c(0.02547394, -0.012162592, -0.02547394,    1,
 			0.02611532, -0.011367422, -0.02611532,    1,
			0.02578871, -0.015207250, -0.02578871,    1),
			nrow=3,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

myres<- nelson_estim(group, bonddata, maturity_spectrum, 
   method, fit, weights, startparam=b,control)

i = 2

plot(myres$yields[[i]][,1],myres$yields[[i]][,2])
kurve <- cbind(myres$yields[[i]][,1],spotrates(method,myres$opt_result[[i]]$par,myres$yields[[i]][,1]))
reihung = order(kurve[,1])
kurve <- kurve[reihung,]
lines(kurve[,1],kurve[,2])

#print(myres)
#summary(myres)
#plot(myres)
