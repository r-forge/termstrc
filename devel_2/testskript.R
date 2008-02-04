#rm(list=ls())

#load("eurobonds.RData")
load("govbonds.RData")

source("methods.R")
source("splines.R")
source("tools.R")
source("nelson.R")

#cubic splines example
#group <- c("GERMANY", "AUSTRIA", "ITALY")
#bonddata <- eurobonds
#matrange <- "all"



#print(x)
#summary(x)
#plot(x)

# nelsonsiegel example
group <- c("GERMANY", "AUSTRIA","FRANCE")
bonddata <- govbonds
matrange <- "all"
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "duration"
control <- list(eval.max=100000, iter.max=500)

b <- matrix(c(0,0,0, 1,
 			0,0,0, 1,
			0,0,0, 1),
			nrow=3,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

#x <- nelson_estim(group, bonddata, matrange, 
#                  method, fit, weights, startparam=b,control)
                 # by <- rbind(x$opt_result$GERMANY$par,x$opt_result$AUSTRIA$par, x$opt_result$ITALY$par)
#y <- nelson_estim(group, bonddata, matrange, method, fit="yields", weights, startparam=by, control)   
#z <- nelson_estim(group, bonddata, matrange, method, fit="prices", weights="duration", startparam=b, control)   

#w <- nelson_estim(group=c("GERMANY"), bonddata <- rm_bond(bonddata,c("DE0001135226",
#"DE0001135275"), gr="GERMANY"), matrange, method, fit="prices", weights="duration", startparam=b, control) 


s <- splines_estim(group, bonddata, matrange)

