#rm(list=ls())

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
#y <- splines_estim(group, bonddata, matrange)

# appliction of bond removal function 
#ISIN <- c("IT0003844534",  "IT0003242747") 
#gr <- rep("ITALY",2)

#testdata <- rm_bond(bonddata,ISIN,gr)
#z <- nelson_estim(group= c("ITALY"), bonddata=testdata, matrange, 
#                  method, fit, weights, startparam=b,control)
#summary(z)
#plot(z,ctype="none", error="price")

#testdata <- rm_bond(bonddata,ISIN=c("DE0001134468"), gr=c("GERMANY"))

#dneu <- nelson_estim(group= c("GERMANY"), bonddata=testdata, matrange,  method, fit, weights, startparam=b,control)