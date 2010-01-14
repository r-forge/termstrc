print.termstrc_ns <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for estimation with spot rate function of:\n")
  cat(switch(x$method,"dl"="Diebold/Li","ns"="Nelson/Siegel","sv"="Svensson"),"\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  parameters <- mapply(function(i) x$opt_result[[i]]$par,seq_along(x$opt_result))
  colnames(parameters) <- names(x$opt_result)
  n_par <- as.character(nrow(parameters))
  rownames(parameters) <- switch(n_par,
          "3"=c("beta_0","beta_1","beta_2"),
          "4"=c("beta_0","beta_1","beta_2","tau_1"),
          "6"=c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")) 
  print.default(parameters)
  cat("\n")
  #x
  }

plot.termstrc_ns <-
  function(x,matrange=c(min(mapply(function(i) min(x$y[[i]][,1]),seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]),seq(x$n_group))))
                        ,multiple=FALSE, expoints=unlist(x$expoints), ctype="spot",
                         errors="none",
                        lwd=2,lty=1,type="l",inset=c(0.8,0.1),ask=TRUE,
                        ...) {
     
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
    
               				
    
    
    cdata <- switch(ctype, "spot" = x$spot,
    					   "forward" = x$forward,
    					   "discount" = x$discount
    					    )
    					   
    cname <- switch(ctype, "spot" = "Zero-coupon yield curve",
    					   "forward" = "Forward rate curve",
    					   "discount" = "Discount factor curve" )
    
    
    # plot all interst rate curves together
    if (multiple) {
    
    plot(x=cdata,multiple=multiple, expoints=expoints,lwd=lwd,type=type,...) }
  
    if (!multiple && ctype %in% c("spot", "forward", "discount")){
        old.par <- par(no.readonly = TRUE)

        if(x$n_group != 1) par(ask=ask)
        
    # plot each interest rate curve seperately
    for (k in seq(x$n_group)  ) 
    	{
    	
    	plot.ir_curve(cdata[[k]], ylim=c(0, max(cdata[[k]][,2]) + 0.01 )*100,
    	xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
             min(ceiling(max(x$y[[k]][,1])),matrange[2])), lwd=lwd,type=type, ...
    	)
    	 
    	title(names(x$opt_result)[k])
    	 
    	if(ctype=="spot") {points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
    		legend("bottom",legend=c("Zero-coupon yield curve","Yield-to-maturity"),
                col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))}
        else 	legend("bottom",legend=cname	,col=c("steelblue"), lty = lty , pch=(-1))

    	
    	}     
 	    on.exit(par(old.par))
     }
    
     # plot spread curves 
    if(ctype == "spread") {plot(x$spread,expoints=expoints,
    	xlim= c(max(floor(samplemat[1]),matrange[1]),
  	min(ceiling(samplemat[2]),matrange[2])),lwd=lwd,
  						    ...)
  	}
    # plot errors 
    if(errors %in% c("price", "yield")){
    	
    	edata <- switch(errors,"price" = x$perrors, "yield"= x$yerrors )
    	if(x$n_group == 1) ask= FALSE
       	for(k in seq(x$n_group)){
    		
     		plot.error(edata[[k]],ask=ask,main=x$group[k],
                           ylab=paste("Error ",paste(errors,"s)",sep=""),sep=" ("),...)
    		
    		legend("bottomright", legend=c(paste("  RMSE",
    		        switch(errors,"price" = round(rmse(x$p[[k]],x$phat[[k]]),4),
                       "yield" = round(rmse(x$y[[k]][,2],x$yhat[[k]][,2]),4)) ,sep=": "),
                        paste("AABSE",switch(errors,"price" = round(aabse(x$p[[k]],x$phat[[k]]),4),
                        "yield" = round(aabse(x$y[[k]][,2],x$yhat[[k]][,2]),4)),sep=": ")),bty="n", inset=inset) 
    		
    	}
    	
      }
    				
   
   
}  

summary.termstrc_ns <- function(object,...) {
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
    sumry <- list(gof,convergencegroup,convergence,startparam=x$startparam)
    names(sumry) <- c("gof","convergencegroup","convergence","startparam")
    class(sumry) <- "summary.termstrc_ns"
    sumry
}

print.summary.termstrc_ns<-
    function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$gof)
    cat("\n")
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Startparameters:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$startparam)
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
