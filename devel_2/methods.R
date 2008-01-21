 
###################################################################
#                print-method for nelson                          #
###################################################################

print.nelson <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for Nelson/Siegel, Svensson estimation:\n")
  cat("\n")
  cat("Method:",x$method,"\n")
  cat("Fitted:",x$fit,"\n")
  cat("Weights:",x$weights,"\n")
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  parameters <- mapply(function(i) x$opt_result[[i]]$par,seq_along(x$opt_result))
  colnames(parameters) <- names(x$opt_result)
  n_par <- as.character(nrow(parameters))
  rownames(parameters) <- switch(n_par,
          "4"=c("beta_0","beta_1","beta_2","tau_1"),
          "6"=c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")) 
  print.default(parameters)
  cat("\n")
  x
  }

###################################################################
#                    plot-method for nelson                       #
###################################################################

plot.nelson <-
  function(x,matrange=c(min(mapply(function(i) min(x$y[[i]][,1]),seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]),seq(x$n_group))))
                        ,pdf=FALSE,multiple=FALSE, expoints=unlist(x$expoints), ctype="spot",
                        lwd=2,
                        ...) {
   
    
     if(pdf) pdf( file="termstrc_results.pdf",... ) else par(ask=TRUE)  
     
     # min and max maturity of all bonds in the sample 
     samplemat <- c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                    max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))) 
  
     # check plot maturity conformity
    if(x$matrange[1] != "all") {
    if(matrange[2]>  x$matrange[2]) { matrange[2] <-  x$matrange[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
   
    if(matrange[1] <  x$matrange[1]) { matrange[1] <-  x$matrange[1]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
    }
   
    if( matrange[2] > samplemat[2]) {matrange[2] <-  samplemat[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
     
    if( matrange[1] < samplemat[1]) {matrange[1] <- samplemat[1]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
    
               				
    
    
    cdata <- switch(ctype, "spot" = x$zcy_curves,
    					   "forward" = x$fwr_curves,
    					   "discount" = x$df_curves )
    					   
    cname <- switch(ctype, "spot" = "Zero-coupon yield curve",
    					   "forward" = "Forward rate curve",
    					   "discount" = "Discount factor curve" )
    
    
    # plot all interst rate curves together
    if (!is.character(x$scurves) && multiple) {
    
    plot(x=cdata,multiple=multiple, expoints=expoints,...) }
 
    
    else 
    {
    
    # plot each interest rate curve seperately
    for (k in seq(x$n_group)  ) 
    	{
    	
    	plot.ir_curve(cdata[[k]], ylim=c(0, max(cdata[[k]][,2]) + 0.01 )*100,
    	xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
             min(ceiling(max(x$y[[k]][,1])),matrange[2]))
    	)
    	 
    	title(names(x$opt_result)[k])
    	 
    	if(ctype=="spot") {points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
    		legend("bottom",legend=c("Zero-coupon yield curve","Yield to maturity"),  			col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))		} else 	legend("bottom",legend=cname	,col=c("steelblue"), lty = 1 , pch=(-1))

    	
    	}     
 	}
     
  # plot spread curves    
  # if (!is.character(x$scurves) && ctype=="spot") {
  # plot(0,0, type="n",
  #  col=(which.max(mapply(function(i) max(x$y[[i]][,1]),
  #                                       seq(x$n_group))[-1]) + 1),lty=1,lwd=lwd,
  #  xlab="Maturity (years)",
  #  ylab="Spread (basis points)",
  #  xlim= c(max(floor(samplemat[1]),matrange[1]),
  #          min(ceiling(samplemat[2]),matrange[2])),
  #  ylim=c(min(x$scurves[,seq(x$n_group)-1]),max(x$scurves[,seq(x$n_group)-1]))*10000)
    
  #  for(k in c(2:x$n_group))
   # { 
   #  lines(x$zcy_curves[1:expoints[k] ,1],x$scurves[1:expoints[k] ,k-1]*10000,col=k,lwd=lwd)
   #  lines(x$zcy_curves[((expoints[k]+1) : nrow(x$zcy_curves) ) ,1],
   #  x$scurves[((expoints+1) : nrow(x$zcy_curves) ) ,k-1]*10000,col=k,lty=5,lwd=lwd)
 #	  } 
   # title("Spread curves")
  # legend("topleft",legend=names(x$opt_result[-1]),col=2:x$n_group,lty=1,lwd=lwd)
   
  # }                    
  
   
   
   if(pdf) dev.off() else par(ask=FALSE) 
   
   
}  

###################################################################
#                 summary-method for nelson                       #
###################################################################

summary.nelson <-
    function(object,...) {
    x <- object
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields","AABSE-Yields")
    convergencegroup <- as.matrix(apply(as.matrix(mapply(function(i) x$opt_result[[i]]$convergence,
                              seq_along(x$opt_result))),1,
                              function(x) if(x==1) "no convergence" else "converged"))
    colnames(convergencegroup) <- "Convergence ()"
    rownames(convergencegroup) <- x$group
    convergence <- as.matrix(mapply(function(i) x$opt_result[[i]]$message,seq_along(x$opt_result)))
    colnames(convergence) <- "Solver message"
    rownames(convergence) <- x$group
    sumry <- list(gof,convergencegroup,convergence)
    names(sumry) <- c("gof","convergencegroup","convergence")
    class(sumry) <- "summary.nelson"
    sumry
}

###################################################################
#                 print-method for summary.nelson                 #
###################################################################

print.summary.nelson <-
    function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$gof)
    cat("\n")
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Convergence information:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$convergencegroup)
    cat("\n")
    print.default(x$convergence)
    cat("\n")
    cat("\n")
    x
}
 
###################################################################
#                print-method for cubic splines                   #
###################################################################

print.cubicsplines <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for Cubic splines estimation:\n")
  cat("\n")
  for(i in seq(x$n_group)) {
  print.default(paste(names(x$alpha)[[i]],":",sep=""))
  names(x$alpha[[i]]) <- paste("alpha",c(seq_along(x$alpha[[i]])))
  print.default(x$alpha[[i]])
  cat("\n")
  x
  }
 }
 
###################################################################
#            summary-method for cubic splines                      #
###################################################################

summary.cubicsplines <-
    function(object,...) {
    x <- object
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields","AABSE-Yields")
    regsumry <- lapply(x$regout,summary)
    for (i in seq(x$n_group)) rownames(regsumry[[i]]$coefficients) <- 
	paste("alpha",c(seq_along(x$alpha[[i]])))
    sumry <- list(gof,regsumry)
    names(sumry) <- c("gof", "regsumry")
    class(sumry) <- "summary.cubicsplines"
    sumry
} 

###################################################################
#            print-method for summary.cubicsplines                #
###################################################################

print.summary.cubicsplines <-
    function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$gof)
    cat("\n")
    x$gof
    
    cat("---------------------------------------------------\n")
    cat("Summary statistics for the fitted models:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$regsumry)
    cat("\n")
    x$regsumry

    
}

###################################################################
#                plot-method for cubic splines                    #
###################################################################

plot.cubicsplines <-
  function(x,matrange =c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group))))
                        ,pdf=FALSE, ...) {
  
     if(pdf) pdf( file="termstrc_results.pdf",... )  else par(ask=TRUE) 
     
     # min and max maturity of all bonds in the sample 
     samplemat <- c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                    max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))) 
  
     # check plot maturity conformity
    if(x$matrange[1] != "all") {
    if(matrange[2]>  x$matrange[2]) { matrange[2] <-  x$matrange[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
    }
   
    if(matrange[1] <  x$matrange[1]) { matrange[1] <-  x$matrange[1]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
    }
   
    if( matrange[2] > samplemat[2]) {matrange[2] <-  samplemat[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
                   				
    # plot each zero cupon yield curve seperately
    for (k in seq(x$n_group)  ) {  
      plot(x$zcy_curves[[k]][,1] ,x$zcy_curves[[k]][,2]*100,
      type="l",
      ylim=c(0,max(x$zcy_curves[[k]][,2],x$zcy_curves[[k]][,3]) + 0.01 )*100,
      xlim=c(max(0,matrange[1]),min(max(x$zcy_curves[[k]][,1]),matrange[2])),
      xlab="Maturity (years)",
      ylab="Percent",
      lwd=2,
      col="steelblue")
      # lower ci         
      lines(x$zcy_curves[[k]][,1],x$zcy_curves[[k]][,3]*100, type="l", lty=3, col="steelblue" )   
      # upper ci 
      lines(x$zcy_curves[[k]][,1],x$zcy_curves[[k]][,4]*100, type="l", lty=3, col="steelblue" )    
      
      title(x$group[k])
      legend("bottomright",legend=c("Zero-coupon yield curve", "95 % Confidence intervall" ,"Yield to 		maturity", "Knot points")
      , col=c("steelblue","steelblue","red", "darkgrey"), lty = c(1,3,-1,2), pch=c(-1,-1,21,-1))
      points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red")
      abline(v=c(x$T[[k]]),lty=2, col="darkgrey") 
      if(pdf == FALSE) par(ask=TRUE) 
    }
    
    
    # plot all zero cupon yield curves together
    if (is.numeric(x$scurves)){
    plot(x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]),
                                 seq(x$n_group)))]][,1],
      		x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]),
                                  seq(x$n_group)))]][,2]*100,
     type="l",
     ylim=c(0,
                max(x$zcy_curves[[which.max(mapply(function(k) 
                max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,2])+ 0.01 )*100,
     xlim=c(max(0,matrange[1]),min(matrange[2],
              max(x$zcy_curves[[which.max(mapply(function(k) 
              max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,1]))),
     xlab="Maturity (years)",
     ylab="Percent",
     lwd=2,
     col=which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]), seq(x$n_group))))
     title("Zero coupon yield curves") 
      
	  for(k in c( (seq(x$n_group))[- which.max(mapply(function(k) 
              max(x$zcy_curves[[k]][,1]), seq(x$n_group))) ]))
	  {lines(x$zcy_curves[[k]][,1] ,
      		x$zcy_curves[[k]][,2]*100, lwd=2,col=k )
		}
	  legend("bottomright",legend=x$group,col=seq(x$n_group), lty=1, lwd=2)	
    }
    
    # plot spread curves    
    if (is.numeric(x$scurves)) {
    matplot(x$zcy_curves[[which.min(mapply(function(k)
             max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,1],
             x$scurves[,1:(x$n_group-1)]*10000, type="l",
    col=2:x$n_group,lty=1,lwd=2,
    xlab="Maturity (years)",
    ylab="Spread (basis points)",
    xlim= c(max(0,matrange[1]),
            min(max(x$zcy_curves[[which.min(mapply(function(k)
            max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,1]),matrange[2])),
    ylim=c(min(x$scurves),max(x$scurves ))*10000)
    title("Spread curves")
    legend("topleft",legend=x$group[-1],col=2:x$n_group,lty=1,lwd=2)
    }                    
  
   if(pdf) dev.off() else par(ask=FALSE) 
}  

###################################################################
#                    plot-method for ir_curve                     #
###################################################################
plot.ir_curve <- function(x,ylim=c(),xlim=c(),lwd=2, type="l",
				xlab="Maturity (years)",ylab="Percent", 
				col="steelblue", ...) {
	
         
      plot(x[,1] ,x[,2]*100,
      type=type,
      ylim=ylim,
      xlim=xlim,
      xlab=xlab,
      ylab=ylab,
      lwd=lwd,
      col=col)
      
}

###################################################################
#                    plot-method for spot_curves                  #
###################################################################

plot.spot_curves <- function(x,multiple= FALSE,
					ylim= c(range(mapply(function(i) range(x[[i]][,2]),seq(x))))*100,xlim=c(),type="l", lty=1, lwd=2, expoints=NULL, ylab= "Zero-coupon yields (%)",
					xlab= "Maturity (years)",main="Zero-coupon yield curves",	 
					
					
					...) {
	
	par(ask=TRUE)
		
	if(multiple) 
	{
		
     plot(x[[which.max(mapply(function(i) max(x[[i]][,1]), seq(x)))]][,1], x[[which.max(mapply(function(i) max(x[[i]][,1]), seq(x)))]][,2]*100, 
      type=type,col=which.max(mapply(function(i) max(x[[i]][,1]), seq(x))),
      lty=lty,lwd=lwd,xlab=xlab,
      ylab=ylab,
      ylim= ylim, ...
      
     )
   
	  for(k in c((seq(x))[-which.max(mapply(function(i) max(x[[i]][,1]), seq(x)))]))
	  
	 {
     lines(x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) else seq(nrow(x[[k]]))),1],x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) else seq(nrow(x[[k]]))),2]*100,col=k,lwd=lwd, ... )
      
      if(is.numeric(expoints))
      {
      lines(x[[k]][((expoints[k]+1):nrow(x[[k]])) ,1],x[[k]][((expoints[k]+1):nrow(x[[k]])),2]*100,col=k,lwd=lwd,lty=5, ... )
      }
           
 	 }
    title(main)
    legend("bottom",legend=names(x),col=seq(x),lty=lty,lwd=lwd)
   }

	else
	{
		for(k in seq(x)) 
		{ 
		plot.ir_curve(x[[k]],...) 
		title(names(x)[k])
      	legend("bottom",legend=main,col=c("steelblue"), lty = 1 , pch=c(-1))
        }		
	}
	par(ask=FALSE)

}

###################################################################
#                    plot-method for fwr_curves                   #
###################################################################

plot.fwr_curves <- function(x,multiple= FALSE,
					ylim= c(range(mapply(function(i) range(x[[i]][,2]),seq(x))))*100,xlim=c(),type="l", lty=1, lwd=2, expoints=NULL, ylab= "Forward rate (%)",
					xlab= "Maturity (years)",main="Forward rate curves",...) 
		{	
			
		plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main,multiple=multiple,expoints=expoints,lty=lty,lwd=lwd, ... )

		}

###################################################################
#                    plot-method for fwr_curves                   #
###################################################################

plot.df_curves <- function(x,multiple= FALSE,
					ylim= c(range(mapply(function(i) range(x[[i]][,2]),seq(x))))*100,xlim=c(),type="l", lty=1, lwd=2, expoints=NULL, ylab= "Discount factor (%)",
					xlab= "Maturity (years)",main="Discount factor curves",...) 
		{	
			
		plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main,multiple=multiple,expoints=expoints,lty=lty,lwd=lwd, ... )

		}
	
	