data(corpbonds)


group <- c("AAA", "A", "BBB")
bonddata <- corpbonds
matrange <- "all"  

x <- splines_estim(group, bonddata, matrange)

print(x)
summary(x)
plot(x)
