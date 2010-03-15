data(govbonds)
data(govbonds)

## Cubic splines estimation
cs_res <- estim_cs(govbonds,c("FRANCE"),matrange=c(0,30))
print(cs_res)
summary(cs_res)

plot(cs_res)

## Pricing errors per bond
plot(cs_res,errors="price",inset=c(.1,.3))
