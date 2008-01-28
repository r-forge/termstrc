rm(list=ls())

load("eurobonds.RData")

source("methods.R")
source("splines.R")
source("tools.R")
source("nelson.R")

# cubic splines example
#group <- c("GERMANY", "AUSTRIA", "ITALY")
#bonddata <- eurobonds
#matrange <- "all"



#print(x)
#summary(x)
#plot(x)

# nelsonsiegel example
group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- c(0,12)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "none"
control <- list(eval.max=100000, iter.max=500)

b <- matrix(c(0,0,0, 1,
 			0,0,0, 1,
			0,0,0, 1),
			nrow=3,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

x <- nelson_estim(group, bonddata, matrange, 
                  method, fit, weights, startparam=b,control)

#print(x)
#summary(x)
#plot(x)
y <- splines_estim(group, bonddata, matrange)
