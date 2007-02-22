data(corpbonds)

countries <- c("AAA","A+","BBB")

bonddata <- corpbonds
maturity_spectrum <- "all"
method="Nelson/Siegel"
fit = "prices"
weights="duration"
control=list(eval.max=100000)

b= matrix(c(0.02794349, -0.015743896, -0.02794349,1,
    		0.02963663, -0.013716390, -0.02963663,1,
    		0.03010361, -0.016545279,-0.03010361,1),
    		nrow=3,ncol=4,byrow=TRUE)

rownames(b) <- countries
colnames(b)<-c("beta0","beta1","beta2","tau1")

myres<- termstrc_estim(countries, bonddata, maturity_spectrum, 
    method, fit, weights, startparam=b,control)   							
print(myres)
summary(myres)
plot(myres)