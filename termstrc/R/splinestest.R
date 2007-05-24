rm(list = ls())
load("eurobonds.RData")
source("splines.R")
source("tools.R")
source("methods.R")

# euro 01
group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- "all" 


myres<- splines_estim(group, bonddata, matrange)

# euro 03

group <- c("GERMANY") 
bonddata <- eurobonds
matrange <- "all"


# corp 01
load("corpbonds.RData")                                                                         
                                                                                                
group <- c("AAA","AA","AA-")                                                                    
bonddata <-  corpbonds                                                                          
matrange <- "all"                                                                                                                                                                              
myres <- splines_estim(group, bonddata, matrange) 

# corp 02 

group <- c("AAA","A+","BBB")
bonddata <-  corpbonds                                                                          
matrange <- "all"                                                                                                                                                                              
myres <- splines_estim(group, bonddata, matrange)                                                