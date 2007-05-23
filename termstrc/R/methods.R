 
###################################################################
#                print-method for nelson                          #
###################################################################

print.nelson <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for Nelson/Siegel, Svensson estimation:\n")
  cat("\n")
  cat("Method:",method,"\n")
  cat("Fitted:",fit,"\n")
  cat("Weights:",weights,"\n")
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  parameters <- mapply(function(i) x$opt_result[[i]]$par,1:length(x$opt_result))
  colnames(parameters) <- names(x$opt_result)
  n_par <- as.character(nrow(parameters))
  rownames(parameters) <- switch(n_par,
          "4"=c("beta_0","beta_1","beta_2","tau_1"),
          "6"=c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")) 
  print.default(parameters)
  cat("\n")
  }

###################################################################
#                    plot-method for nelson                       #
###################################################################

plot.nelson <-
  function(x,matrange = c( min(mapply(function(i) min(x$y[[i]][,1]), 1:x$n_group)),
                           max(mapply(function(i) max(x$y[[i]][,1]), 1:x$n_group)))
                          ,pdf=FALSE, ...) {
   
     if(pdf== TRUE) pdf( file="termstrc_ouput.pdf",... )
     
     
     # min and max maturity of all bonds in the sample 
     samplemat <- c( min(mapply(function(i) min(x$y[[i]][,1]), 1:x$n_group)),
                           max(mapply(function(i) max(x$y[[i]][,1]), 1:x$n_group))) 
  
     # check plot maturity conformity
    
    if(x$matrange != "all") {
    if(matrange[2]>  x$matrange[2]) { matrange[2] <-  x$matrange[2]
     warning("The maximum plot maturity range exceeds the choosen maximum maturity considered for the estimation")}
   
    if(matrange[1] <  x$matrange[1]) { matrange[1] <-  x$matrange[1]
     warning("The minium plot maturity range undercuts the choosen minimum maturity considered for the estimation")}
    }
   
    if( matrange[2] > samplemat[2]) {matrange[2] <-  samplemat[2]
     warning("The maximum plot maturity range has been set to the maxium maturity of the bond with the longest maturity" )
     }
     
    if( matrange[1] < samplemat[1]) {matrange[1] <- samplemat[1]
     warning("The minium plot maturity range has been set to the minium maturity of the bond with the shortest maturity" )
     }
                       				
    # plot each yield curve seperately
    for (k in 1:x$n_group  ) {
         
      plot(x$zcy_curves[,1] ,x$zcy_curves[,k+1]*100,
      type="l",
      ylim=c(0, max(x$y[[k]][,2]) + 0.01 )*100,
      xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]), min(ceiling(max(x$y[[k]][,1])),matrange[2])),
      xlab="Maturities (in years) ",
      ylab="Percent",
      lwd=2,
      col="steelblue")
      title(names(x$opt_result)[k])
      legend("bottomright",legend=c("Zero-coupon yield curve","Yields"),col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))
      grid()
      points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
      if(pdf == FALSE) par(ask=TRUE) 
    }
     
     
      
    # plot all zero coupon yield curves together
    plot(x$zcy_curves[,1], x$zcy_curves[, which.max(mapply(function(i) max(x$y[[i]][,1]), 1:x$n_group)) +1 ]*100, type="l",
    col=which.max(mapply(function(i) max(x$y[[i]][,1]), 1:x$n_group)),lty=1,lwd=2,
    xlab="Maturities (in years)",
    ylab="Zero-Coupon yields (in percent)",
    xlim=c(max(floor(samplemat[1]),matrange[1]), min(ceiling(samplemat[2]),matrange[2])),
    ylim= c(0,max(x$zcy_curves[,2:x$n_group+1] ))*100)
    
   
	  for(k in c( (1:x$n_group)[- which.max(unlist(lapply(x$y,max)))]))
	  {
     spoint <- which(x$zcy_curves[,1] > unlist(lapply(x$y,max))[k])[1] 
     lines(x$zcy_curves[1:spoint ,1],x$zcy_curves[1:spoint ,k+1]*100,col=k, lwd=2)
     lines(x$zcy_curves[((spoint+1) : nrow(x$zcy_curves) ) ,1],
     x$zcy_curves[((spoint+1) : nrow(x$zcy_curves) ) ,k+1]*100,col=k, lty=5, lwd=2)
 	  } 
    title("Zero-coupon yield curves")
    legend("bottomright",legend=names(x$opt_result),col=1:x$n_group,lty=1,lwd=2)
    grid()
    
    # plot spread curves    
    if ( is.character(x$scurves) == FALSE ) {
    matplot(x$zcy_curves[,1], x$scurves[, which.max(unlist(lapply(x$y,max)))-1]*10000, type="l",
    col=which.max(unlist(lapply(x$y,max))),lty=1,lwd=2,
    xlab="Maturities (in years)",
    ylab="Spread (in bps)",
    xlim= c(max(floor(samplemat[1]),matrange[1]), min(ceiling(samplemat[2]),matrange[2])),
    ylim=c(min(x$scurves[,1:x$n_group-1]),max(x$scurves[,1:x$n_group-1] ))*10000)
    
    for(k in c((1:x$n_group))[c(-1,- which.max(unlist(lapply(x$y,max))))])
	  {
     spoint <- which(x$zcy_curves[,1] > unlist(lapply(x$y,max))[k])[1] 
     
     lines(x$zcy_curves[1:spoint ,1],x$scurves[1:spoint ,k-1]*10000,col=k, lwd=2)
     lines(x$zcy_curves[((spoint+1) : nrow(x$zcy_curves) ) ,1],
     x$scurves[((spoint+1) : nrow(x$zcy_curves) ) ,k-1]*10000,col=k, lty=5, lwd=2)
 	  } 
    title("Spread curves")
    legend("topleft",legend=names(x$opt_result[-1]),col=2:x$n_group,lty=1,lwd=2)
    grid()
    
    
    }                    
  
   if(pdf== TRUE) dev.off()
}  

###################################################################
#                 summary-method for nelson                       #
###################################################################

summary.nelson <-
    function(x) {
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),1:x$n_group)
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),1:x$n_group)
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2],x$yhat[[i]][,2]),1:x$n_group)
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2],x$yhat[[i]][,2]),1:x$n_group)
    
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields","AABSE-Yields")
    cat("---------------------------------------------------\n")
    cat("Goodness of fit tests:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(gof)
    cat("\n")
    cat("\n")
    convergencegroup <- as.matrix(apply(as.matrix(mapply(function(i) x$opt_result[[i]]$convergence,
                              1:length(x$opt_result))),1,
                              function(x) if(x==1) "no convergence" else "converged"))
    colnames(convergencegroup) <- "Convergence ()"
    rownames(convergencegroup) <- x$group
    convergence <- as.matrix(mapply(function(i) x$opt_result[[i]]$message,1:length(x$opt_result)))
    colnames(convergence) <- "Solver message"
    rownames(convergence) <- x$group
    cat("---------------------------------------------------\n")
    cat("Convergence information:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(convergencegroup)
    cat("\n")
    print.default(convergence)
    cat("\n")
    cat("\n")
   
} 
###################################################################
#                print-method for cubic splines                   #
###################################################################

print.cubicsplines <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for Cubic splines estimation:\n")
  cat("\n")
  for(i in 1:x$n_group) {
  print.default(paste(names(x$alpha)[[i]],":",sep=""))
  names(x$alpha[[i]]) <- paste("alpha",c(1:length(x$alpha[[i]])))
  print.default(x$alpha[[i]])
  cat("\n")
  }
  
 
 }
 
 ###################################################################
#            summary-method for cubic splines                      #
###################################################################

summary.cubicsplines <-
    function(x) {
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),1:x$n_group)
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),1:x$n_group)
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2],x$yhat[[i]][,2]),1:x$n_group)
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2],x$yhat[[i]][,2]),1:x$n_group)
    
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields","AABSE-Yields")
    cat("---------------------------------------------------\n")
    cat("Goodness of fit tests:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(gof)
    cat("\n")
   
} 

###################################################################
#                plot-method for cubic splines                    #
###################################################################

plot.cubicsplines <-
  function(x,matrange = c( min(mapply(function(i) min(x$y[[i]][,1]), 1:x$n_group)),
                           max(mapply(function(i) max(x$y[[i]][,1]), 1:x$n_group)))
                          ,pdf=FALSE, ...) {
   
     if(pdf== TRUE) pdf( file="termstrc_ouput.pdf",... )
     
     
     # min and max maturity of all bonds in the sample 
     samplemat <- c( min(mapply(function(i) min(x$y[[i]][,1]), 1:x$n_group)),
                           max(mapply(function(i) max(x$y[[i]][,1]), 1:x$n_group))) 
  
     # check plot maturity conformity
    
    if(x$matrange != "all") {
    if(matrange[2]>  x$matrange[2]) { matrange[2] <-  x$matrange[2]
     warning("The maximum plot maturity range exceeds the choosen maximum maturity considered for the estimation")}
   
    if(matrange[1] <  x$matrange[1]) { matrange[1] <-  x$matrange[1]
     warning("The minium plot maturity range undercuts the choosen minimum maturity considered for the estimation")}
    }
   
    if( matrange[2] > samplemat[2]) {matrange[2] <-  samplemat[2]
     warning("The maximum plot maturity range has been set to the maxium maturity of the bond with the longest maturity" )
     }
     
    if( matrange[1] < samplemat[1]) {matrange[1] <- samplemat[1]
     warning("The minium plot maturity range has been set to the minium maturity of the bond with the shortest maturity" )
     }
                       				
    # plot each zero cupon yield curve seperately
    for (k in 1:x$n_group  ) {  
      plot(x$zcy_curves[[k]][,1] ,x$zcy_curves[[k]][,2]*100,
      type="l",
      ylim=c(0, (max(x$zcy_curves[[k]][,2]) + 0.01 )*100),
      xlim=c(0,max(x$zcy_curves[[k]][,1])),
      xlab="Maturities (in years)",
      ylab="Percent",
      lwd=2,
      col="steelblue")
      title(names(x$n_group)[k])
      legend("bottomright",legend=c("Zero-coupon yield curve","Yields"),col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))
      grid()
      points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red")
      abline(v=c(x$T[[k]]),lty=2, col="grey") 
      par(ask=TRUE) 
    }
    
    
    # plot all zero cupon yield curves together
    
      plot(x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]), 1:x$n_group))]][,1] ,
      		x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]), 1:x$n_group))]][,2]*100,
      type="l",
      ylim=c(0, (max(x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]), 1:x$n_group))]][,2])+ 0.01 )*100),
      xlim=c(0,max(x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]), 1:x$n_group))]][,1])),
      xlab="Maturities (in years)",
      ylab="Percent",
      lwd=2,
      col=which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]), 1:x$n_group)))
      
      title("Zero coupon yield curves") 
      
	  for(k in c( (1:x$n_group)[- which.max(unlist(lapply(x$y,max)))]))
	  {lines(x$zcy_curves[[k]][,1] ,
      		x$zcy_curves[[k]][,2]*100, lwd=2,col=k )
		}
		
	  legend("bottomright",legend=x$group,col=1:x$n_group, lty=1, lwd=2)	  #grid()
      
      
    

    
     
        
    # plot spread curves    
    #if ( is.character(x$scurves) == FALSE ) {
    #matplot(x$zcy_curves[,1], x$scurves[, which.max(unlist(lapply(x$y,max)))-1], type="l",
    #col=which.max(unlist(lapply(x$y,max))),lty=1,lwd=2,
    #xlab="Maturities",
    #ylab="Spread",
    #xlim= c(max(floor(samplemat[1]),matrange[1]), min(ceiling(samplemat[2]),matrange[2])),
    #ylim=c(min(x$scurves[,1:x$n_group-1]),max(x$scurves[,1:x$n_group-1] )))
    
    #for(k in c((1:x$n_group))[c(-1,- which.max(unlist(lapply(x$y,max))))])
	#  {
    # spoint <- which(x$zcy_curves[,1] > unlist(lapply(x$y,max))[k])[1] 
     
    # lines(x$zcy_curves[1:spoint ,1],x$scurves[1:spoint ,k-1],col=k, lwd=2)
    # lines(x$zcy_curves[((spoint+1) : nrow(x$zcy_curves) ) ,1],
    # x$scurves[((spoint+1) : nrow(x$zcy_curves) ) ,k-1],col=k, lty=5, lwd=2)
 	#  } 
    #title("Spread curves")
    #legend("topleft",legend=names(x$alpha[-1]),col=2:x$n_group,lty=1,lwd=2)
    #grid()
    
    
    #}                    
  
   if(pdf== TRUE) dev.off()
}  
