
###################################################################
#               Dynamic term structure estimation                 #
###################################################################



###################################################################
#               Parameters extractor function                     #
###################################################################



###################################################################
#          summary-method for dyntermstrc parameters              #
###################################################################



###################################################################
#              print-method for summary.param                     #
###################################################################




###################################################################
#              summary-method for dyntermstrc                     #
###################################################################




###################################################################
#              print-method for summary.dyntermstrc               #
###################################################################






###################################################################
#                plot-method for dyntermstrc                      #
###################################################################



###################################################################
#              Pricing errors extractor function                  #
###################################################################

perrors <- function(x) {
  # pricing errors
  errors <- t(mapply(function(i) x[[i]]$perrors[[1]][,2], seq(length(x))))
  # maturities at start date (needed for spacing plot method)
  ttmat <- x[[length(x)]]$y[[1]][,1]         
  perrors <- list(errors=errors,ttmat=ttmat)   
  class(perrors) <- "3derror"
  perrors
}


###################################################################
#                plot-method for 3derror                          #
###################################################################

# TODO: switch for yield errors

plot.3derror <- function(x, ... ) {
  old.par <- par(no.readonly = TRUE) 

  # errors <- switch(type,
  #                  "price" = perrors(x),
  #                  "yield" = yerrors(x))

  if (is.list(x)) {             # estimation error as list
    persp(seq(length(x$ttmat)), # use "x$ttmat" for maturity spacing => problems if bonds with same maturity
          seq(nrow(x$errors)),
          t(x$errors),
          theta = -35,
          phi = 30, expand = 0.6, col = "lightgreen",
          ltheta = 120, shade = 0.55, ticktype = "detailed",
          xlab="Maturity",zlab= " Price error",ylab="Time",
          box=TRUE,border=NA)
  } else                        # forecast error as matrix        
    persp(seq(ncol(x)),seq(nrow(x)),t(x), theta = -35, phi = 30,
          expand = 0.6, col = "lightgreen", ltheta = 120,
          shade = 0.55, ticktype = "detailed", xlab="Bond",
          zlab= " Price error",ylab="Time",box=TRUE,border=NA)

  on.exit(par(old.par))  
}


 
###################################################################
#              Yields extractor function                          #
###################################################################

yields <- function(x) {
  yields <- t(mapply(function(i) x[[i]]$y[[1]][,2], seq(length(x))))                  
  yields
}









