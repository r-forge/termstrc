data(eurobonds)

countries <- c("GERMANY")
bonddata <- eurobonds
maturity_spectrum <- "all"
method = "Nelson/Siegel"
fit = "prices"
weights = "duration"
control=list(eval.max=100000)

b<-matrix(c(0.02547394, -0.012162592, -0.02547394, 1),nrow=1,ncol=4,byrow=TRUE)
			
rownames(b)<-countries

colnames(b)<-c("beta0","beta1","beta2","tau1")

myres<- termstrc_estim(countries, bonddata, maturity_spectrum, 
    method, fit, weights, startparam=b,control)
    
print(myres)
summary(myres)
plot(myres,spread_curves=FALSE)
