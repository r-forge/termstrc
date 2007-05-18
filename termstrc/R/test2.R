rm(list = ls())
load("eurobonds.RData")
source("nelson.R")
source("tools.R")

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
maturity_spectrum <- c(2,15)
method = "Svensson"
fit = "yields"
weights = "duration"
control=list(eval.max=100000)

b<-matrix(c(0.04294523, -0.01149559, -0.03321905,  1.136986, -0.01252171,  9.643836,
			0.03848141, -0.01378969, -0.01234509,  2.167123,  0.00566015,  8.668493,
			0.00000000,  0.02159121,  0.10187099, 28.728767,  0.01342987,  2.213699),
			nrow=3,ncol=6,byrow=TRUE)

rownames(b)<-group
colnames(b)<-c("beta0","beta1","beta2","tau1","beta3","tau2")

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
