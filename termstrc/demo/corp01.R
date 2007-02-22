data(corpbonds)

countries = c("AAA","AA","AA-")
bonddata = corpbonds
maturity_spectrum = c(4,15)
method = "Svensson"
fit = "prices"
weights = "duration"
control=list(eval.max=100000)

b=matrix(c(1.10966150, -1.08396184, -0.35424537,  7.958904, -2.64900449, 29.863014,
 			0.07198431, -0.04077098, -0.04122995,  2.202740, -0.07079683, 11.030137,
 			0.26592507, -0.24258357, -0.08361875,  4.835616, -0.54029496, 18.964384),
			nrow=3,ncol=6,byrow=TRUE)
			
rownames(b)<-countries
colnames(b)<-c("beta0","beta1","beta2","tau1","beta3","tau2")			

myres<- termstrc_estim(countries, bonddata, maturity_spectrum, 
    method, fit, weights, startparam=b,control)
    
print(myres)
summary(myres)
plot(myres,4,15)


