# euro 01
rm(list = ls())
load("eurobonds.RData")
source("nelson.R")
source("tools.R")
source("methods.R")

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- c(2,10)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "none"
control <- list(eval.max=100000)

b <- matrix(c(0.02547394, -0.012162592, -0.02547394,    1,
 			0.02611532, -0.011367422, -0.02611532,    1,
			0.02578871, -0.015207250, -0.02578871,    1),
			nrow=3,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

myres <- nelson_estim(group, bonddata, matrange, 
   method, fit, weights, startparam=b,control)

# euro 02

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- "all"
method <- "Svensson"
fit <- "prices"
weights <- "duration"
control <- list(eval.max=100000)

b <- matrix(c(0.04294523, -0.01149559, -0.03321905,  1.136986, -0.01252171,  9.643836,
			0.03848141, -0.01378969, -0.01234509,  2.167123,  0.00566015,  8.668493,
			0.00000000,  0.02159121,  0.10187099, 28.728767,  0.01342987,  2.213699),
			nrow=3,ncol=6,byrow=TRUE)

rownames(b) <- group
colnames(b) <- c("beta0","beta1","beta2","tau1","beta3","tau2")

myres <- nelson_estim(group,bonddata, matrange, 
    method, fit, weights, startparam=b,control)	
    
#euro 03

group <- c("GERMANY")
bonddata <- eurobonds
matrange <- c(2,10)
method <- "Nelson/Siegel"
fit <- "yields"
weights <- "duration"
control <- list(eval.max=100000)

b <- matrix(c(0.02547394, -0.012162592, -0.02547394, 1),nrow=1,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

myres <- nelson_estim(group, bonddata, matrange, 
    method, fit, weights, startparam=b,control)
    
# corp 01 

load("corpbonds.RData")

group <- c("AAA","AA","AA-")
bonddata <-  corpbonds
matrange <- "all"
method <- "Svensson"
fit <- "yields"
weights <- "duration"
control <- list(eval.max=100000)

b <- matrix(c(1.10966150, -1.08396184, -0.35424537,  7.958904, -2.64900449, 29.863014,
 			0.07198431, -0.04077098, -0.04122995,  2.202740, -0.07079683, 11.030137,
 			0.26592507, -0.24258357, -0.08361875,  4.835616, -0.54029496, 18.964384),
			nrow=3,ncol=6,byrow=TRUE)
			
rownames(b) <- group
colnames(b) <- c("beta0","beta1","beta2","tau1","beta3","tau2")			

myres <- nelson_estim(group, bonddata, matrange, 
    method, fit, weights, startparam=b,control)

#corp 02 
group <- c("AAA","A+","BBB")

bonddata <- corpbonds
matrange <- c(2,10)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "duration"
control <- list(eval.max=100000)

b <- matrix(c(0.02794349, -0.015743896, -0.02794349,1,
    		0.02963663, -0.013716390, -0.02963663,1,
    		0.03010361, -0.016545279,-0.03010361,1),
    		nrow=3,ncol=4,byrow=TRUE)

rownames(b) <- group
colnames(b)<- c("beta0","beta1","beta2","tau1")

myres <- nelson_estim(group, bonddata, matrange, 
    method, fit, weights, startparam=b,control)   
    
    
group <- c("AAA")                                                                  		    
bonddata <-  corpbonds                                                                    
matrange <- c(2,10)                                                                         
method <- "Svensson"                                                                      
fit <- "yields"                                                                           
weights <- "duration"                                                                     
control <- list(eval.max=100000)                                                          
                                                                                          
b <- matrix(c(1.10966150, -1.08396184, -0.35424537,  7.958904, -2.64900449, 29.863014     
 			),nrow=1,ncol=6,byrow=TRUE)                                                         
			                                                                                    
rownames(b) <- group                                                                      
colnames(b) <- c("beta0","beta1","beta2","tau1","beta3","tau2")			                      
                                                                                          
myres <- nelson_estim(group, bonddata, matrange,                                          
    method, fit, weights, startparam=b,control)                                           