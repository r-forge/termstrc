data(eurobonds)

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- c(2,12)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "none"
control <- list(eval.max=100000)

b <- matrix(c(0,0,0, 1,
 			0,0,0, 1,
			0,0,0, 1),
			nrow=3,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

x <- nelson_estim(group, bonddata, matrange, 
                  method, fit, weights, startparam=b,control)

print(x)
summary(x)
plot(x)

