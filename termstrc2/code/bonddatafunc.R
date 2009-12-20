
## Limit a bond data set to a certain maturity range 
maturity_range <- 
  function(bonddata,lower,upper) {

m <- lapply(bonddata,create_maturities_matrix)
colmax <- function(m) apply(m,2,max)
m_max <- lapply(m,colmax)

bonds_in_range <- function(group) names(group[which(
                group>lower & group<upper)])

# list with ISINs of bonds in range
isins_range <- lapply(m_max,bonds_in_range)
index_set <- which(unlist(lapply(isins_range,length))>0)

# list with positions of bonds in isins_range
isins_range_pos <- list()

for (i in seq_along(bonddata)) {
        isins_range_pos[[i]] <- which(bonddata[[i]]$ISIN %in% 
            isins_range[[i]])
    }
names(isins_range_pos) <- names(bonddata)

# first part of bonddata for filtering

N <- which(names(bonddata[[1]]) %in% c("ISIN","MATURITYDATE","STARTDATE","COUPONRATE","PRICE","ACCRUED"))

first <- function(lst) lst[N]
filtered <- lapply(bonddata,first)

bonddata_range <- list()
for (i in seq_along(bonddata)) {
bonddata_range[[i]] <- as.list(as.data.frame(filtered[[i]])
                       [isins_range_pos[[i]],])
# convert to character                       
bonddata_range[[i]][["ISIN"]] <- as.character(bonddata_range[[i]][["ISIN"]])
bonddata_range[[i]][["MATURITYDATE"]] <- as.character(bonddata_range[[i]][["MATURITYDATE"]])
bonddata_range[[i]][["STARTDATE"]] <- as.character(bonddata_range[[i]][["STARTDATE"]])
}
names(bonddata_range) <- names(bonddata)

# list with positions of cashflows in isins_range

isins_range_pos <- list()
for (i in seq_along(bonddata)) {
isins_range_pos[[i]] <- which(bonddata[[i]][["CASHFLOWS"]]
                        [["ISIN"]]%in%isins_range[[i]])
}
names(isins_range_pos) <- names(bonddata)

for (i in seq_along(bonddata)) {
CASHFLOWS <- as.list(as.data.frame(bonddata[[i]]
             [["CASHFLOWS"]])[isins_range_pos[[i]],])
CASHFLOWS$ISIN <- as.character(CASHFLOWS$ISIN)
CASHFLOWS$DATE <- as.character(CASHFLOWS$DATE)             
bonddata_range[[i]][["CASHFLOWS"]] <- list()
bonddata_range[[i]][["CASHFLOWS"]] <- CASHFLOWS
}
names(bonddata_range) <- names(bonddata)
bonddata_range

# add TODAY from bonddata
for (i in seq_along(bonddata)) {
bonddata_range[[i]][["TODAY"]] <- bonddata[[i]][["TODAY"]]
}

# delete countries where no bonds are available
bonddata_range <- bonddata_range[index_set]
bonddata_range
}




  
      


###################################################################
#                    Bond removal function                        #
###################################################################

rm_bond <- function(bdata,ISIN,gr){
    cf_isin_index <- which(bdata[[gr]]$CASHFLOWS$ISIN %in% ISIN)
 	isin_index <- which(bdata[[gr]]$ISIN %in% ISIN)	

    	bdata[[gr]]$ISIN <-  bdata[[gr]]$ISIN[-isin_index]
    	bdata[[gr]]$MATURITYDATE <- bdata[[gr]]$MATURITYDATE[-isin_index]
    	bdata[[gr]]$STARTDATE <- bdata[[gr]]$STARTDATE[-isin_index]
    	bdata[[gr]]$COUPONRATE <- bdata[[gr]]$COUPONRATE[-isin_index]
    	bdata[[gr]]$PRICE <- bdata[[gr]]$PRICE[-isin_index]
    	bdata[[gr]]$ACCRUED <- bdata[[gr]]$ACCRUED[-isin_index]

		bdata[[gr]]$CASHFLOWS$ISIN <- bdata[[gr]]$CASHFLOWS$ISIN[-cf_isin_index]
		bdata[[gr]]$CASHFLOWS$CF <- bdata[[gr]]$CASHFLOWS$CF[-cf_isin_index]
		bdata[[gr]]$CASHFLOWS$DATE <- bdata[[gr]]$CASHFLOWS$DATE[-cf_isin_index]
	
	bdata
}



###################################################################
#           Bonddata preprocessing function                       #
###################################################################
prepro_bond <- function(group,
           bonddata,
           matrange="all",
           bpeq="dirty"){

  # select given group from bonddata
  bonddata <- bonddata[group]

  # select data according to chosen maturity range
  if (length(matrange)==1) {bonddata <- bonddata }else
   {bonddata <- maturity_range(bonddata,matrange[1],matrange[2]) }

  # number of groups 
  n_group <- length(bonddata) 
  
  # group sequence
  sgroup <- seq(n_group)
    
  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)

  # create cashflows matrix including dirty price (needed for bond yield calculation)
  cf_pd <- mapply(function(k) create_cashflows_matrix(bonddata[[k]],include_price=TRUE,ai=TRUE),
                 sgroup,SIMPLIFY=FALSE)
 
  # create cashflows matrix including clean price
  if(bpeq=="clean") {
    cf_p <- mapply(function(k) create_cashflows_matrix(bonddata[[k]],include_price=TRUE,ai=FALSE),
                 sgroup,SIMPLIFY=FALSE)
  } else {
    cf_p = cf_pd
  }
  
  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # create maturities matrix including zeros (needed for bond yield calculation)
  m_p <- mapply(function(k) create_maturities_matrix(bonddata[[k]],include_price=TRUE),
                sgroup,SIMPLIFY=FALSE)
  
  # calculate dirty prices
  pd <- mapply(function(k) bonddata[[k]]$PRICE + bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)

  # extract clean prices (new: if-condition)
  if(bpeq=="clean") {
    p <- mapply(function(k) bonddata[[k]]$PRICE,sgroup,SIMPLIFY=FALSE)
  } else {
    p <- pd
  }
  
  # extract accrued interest
  ac <- mapply(function(k) bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)

  # browser()
  # assign ISIN 
  for(k in sgroup) {names(pd[[k]]) <- bonddata[[k]]$ISIN
                    names(p[[k]]) <- bonddata[[k]]$ISIN
                    names(ac[[k]]) <- bonddata[[k]]$ISIN
                 
                  }
  
  # index for ordering
  positions <- mapply(function(k) order(apply(m[[k]],2,max)),sgroup,SIMPLIFY=FALSE)
  
  
  # order matrices 
  cf <- mapply(function(k) cf[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  cf_p <- mapply(function(k) cf_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  cf_pd <- mapply(function(k) cf_pd[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m <- mapply(function(k) m[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m_p <- mapply(function(k) m_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  pd <- mapply(function(k) pd[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  p <- mapply(function(k) p[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  ac <- mapply(function(k) ac[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  # browser()
  # calculate bond yields	
  y <- mapply(function(k) bond_yields(cf_pd[[k]],m_p[[k]]),
                   sgroup,SIMPLIFY=FALSE)

  # calculate bond yield based on the clean price

  yc <- mapply(function(k) bond_yields(cf_p[[k]],m_p[[k]]),sgroup,SIMPLIFY=FALSE)
  # calculate duration   
  duration <- mapply(function(k) duration(cf_pd[[k]],m_p[[k]],y[[k]][,2]),
                   sgroup,SIMPLIFY=FALSE)

  durationc <- mapply(function(k) duration(cf_p[[k]],m_p[[k]],yc[[k]][,2]),
                   sgroup,SIMPLIFY=FALSE)

res <- list(n_group=n_group,sgroup=sgroup,positions=positions,cf=cf,cf_p=cf_p,cf_pd=cf_pd,m=m,m_p=m_p,pd=pd,p=p,ac=ac,y=y,yc=yc,duration=duration,durationc=durationc,timestamp=bonddata[[1]]$TODAY)
res
}




###################################################################
#           Bonddata postprocessing function                      #
###################################################################




postpro_bond <- function(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,method){
 

 
 # theoretical bond prices with estimated parameters
 phat <- mapply(function(k) bond_prices(method,opt_result[[k]]$par,
       m[[k]],cf[[k]])$bond_prices,sgroup,SIMPLIFY=FALSE)

 
  # price errors
  perrors <- mapply(function(k) cbind(y[[k]][,1],phat[[k]] - p[[k]]),sgroup,SIMPLIFY=FALSE)     
 
  for (k in sgroup) class(perrors[[k]]) <- "error"

  # calculate estimated dirty prices pdhat = phat + a
  pdhat <- mapply(function(k) phat[[k]]+ac[[k]],sgroup,SIMPLIFY=FALSE)
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-pdhat[[k]],cf[[k]]),m_p[[k]]),sgroup,SIMPLIFY=FALSE)
  
  # yield errors
  yerrors <- mapply(function(k) cbind(y[[k]][,1], yhat[[k]][,2] - y[[k]][,2]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(yerrors[[k]]) <- "error"

  
  # maturity interval
  t <- seq(round(min(mapply(function(i) min(y[[i]][,1]), sgroup)),2),
                           ceiling(max(mapply(function(i) max(y[[i]][,1]), sgroup))),0.01)
  
  
  # calculate zero coupon yield curves  
  zcy_curves <- switch(method,
              "Nelson/Siegel" = mapply(function(k)
		            cbind(t,nelson_siegel(opt_result[[k]]$par,t)),sgroup, SIMPLIFY=FALSE),
              "Svensson" = mapply(function(k) 
              		cbind(t,svensson(opt_result[[k]]$par,t)),sgroup,SIMPLIFY=FALSE),
              "Diebold" = mapply(function(k)
		            cbind(t,diebold(opt_result[[k]]$par,t)),sgroup, SIMPLIFY=FALSE))
              
                       
              		 

  for (k in sgroup) class(zcy_curves[[k]]) <- "ir_curve"
  class(zcy_curves) <- "spot_curves"
                                      
  # calculate spread curves              	 
   if(n_group != 1) {  
   s_curves <- mapply(function(k) cbind(t,(zcy_curves[[k]][,2] - zcy_curves[[1]][,2])),sgroup,
   					SIMPLIFY=FALSE)
    } else s_curves = "none"
   
   for (k in sgroup) class(s_curves[[k]]) <- "ir_curve" 
   class(s_curves) <- "s_curves"
    
  # calculate extrapolation point                        
  expoints <- mapply(function(k) which(zcy_curves[[k]][,1] > 
                 mapply(function(i) max(y[[i]][,1]), seq(n_group))[k])[1],sgroup, SIMPLIFY=FALSE )  
        
  # calculate forward rate curves 
  fwr_curves <- switch(method,
              "Nelson/Siegel" = mapply(function(k)
		            cbind(t,fwr_ns(opt_result[[k]]$par,t)),sgroup, SIMPLIFY=FALSE),
              "Svensson" = mapply(function(k) 
              		cbind(t,fwr_sv(opt_result[[k]]$par,t)),sgroup,SIMPLIFY=FALSE),
               "Diebold" = mapply(function(k)
		            cbind(t,fwr_db(opt_result[[k]]$par,t)),sgroup, SIMPLIFY=FALSE))
                      
                   
                      
  for (k in sgroup) class(fwr_curves[[k]]) <- "ir_curve"
  class(fwr_curves) <- "fwr_curves"

  
  # calculate discount factor curves 
  df_curves <- mapply(function(k) cbind(zcy_curves[[k]][,1],exp(-zcy_curves[[k]][,1]*
                                        zcy_curves[[k]][,2])),sgroup,SIMPLIFY=FALSE)
  
   for (k in sgroup) class(df_curves[[k]]) <- "ir_curve"
   class(df_curves) <- "df_curves"

  res <- list(phat=phat,perrors=perrors,pdhat=pdhat,yhat=yhat,yerrors=yerrors,t=t,zcy_curves=zcy_curves,
              s_curves=s_curves,expoints=expoints,fwr_curves=fwr_curves,df_curves=df_curves,opt_result=opt_result)
  res

}

###################################################################
#           Bonddata diagnostics function                         #
###################################################################

diag_bond <- function(bonddata)  {

  # number of groups 
  n_group <- length(bonddata) 
  
  # group sequence
  sgroup <- seq(n_group) 

  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)

  # create cashflows matrix including dirty price (needed for bond yield calculation)
  cf_pd <- mapply(function(k) create_cashflows_matrix(bonddata[[k]],include_price=TRUE,ai=TRUE),
                 sgroup,SIMPLIFY=FALSE)
  
  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # create maturities matrix including zeros (needed for bond yield calculation)
  m_p <- mapply(function(k) create_maturities_matrix(bonddata[[k]],include_price=TRUE),
                sgroup,SIMPLIFY=FALSE)
  
  # calculate dirty prices
  pd <- mapply(function(k) bonddata[[k]]$PRICE + bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)
  p <- mapply(function(k) bonddata[[k]]$PRICE,sgroup,SIMPLIFY=FALSE)
  # assign ISIN 
  for(k in sgroup) names(pd[[k]]) <- bonddata[[k]]$ISIN
  for(k in sgroup) names(p[[k]]) <- bonddata[[k]]$ISIN
  # index for ordering
  positions <- mapply(function(k) order(apply(m[[k]],2,max)),sgroup,SIMPLIFY=FALSE)
  
  
  # order matrices 
  cf <- mapply(function(k) cf[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  cf_pd <- mapply(function(k) cf_pd[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m <- mapply(function(k) m[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m_p <- mapply(function(k) m_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  pd <- mapply(function(k) pd[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  p <- mapply(function(k) p[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  
  # calculate bond yields	
  y <- mapply(function(k) bond_yields(cf_pd[[k]],m_p[[k]]),
              sgroup,SIMPLIFY=FALSE)
  res <- list(cf_pd=cf_pd,cf=cf,m_p=m_p,m=m,y=y,p=p,pd=pd)
  for ( i in 1:length(res)) names(res[[i]]) <- names(bonddata)
  res
  }

###################################################################
#               Remove bonds from dynamic data set                #
###################################################################

dyn_rm_bond <- function(bdata, ISIN) {
  for (i in seq(length(bdata))) {
     cf_isin_index <- which(bdata[[i]]$CASHFLOWS$ISIN %in% ISIN)
     isin_index <- which(bdata[[i]]$ISIN %in% ISIN)
     bdata[[i]]$ISIN <- bdata[[i]]$ISIN[-isin_index]
     bdata[[i]]$MATURITYDATE <- bdata[[i]]$MATURITYDATE[-isin_index]
     bdata[[i]]$STARTDATE <- bdata[[i]]$STARTDATE[-isin_index]
     bdata[[i]]$COUPONRATE <- bdata[[i]]$COUPONRATE[-isin_index]
     bdata[[i]]$PRICE <- bdata[[i]]$PRICE[-isin_index]
     bdata[[i]]$ACCRUED <- bdata[[i]]$ACCRUED[-isin_index]
     bdata[[i]]$RY <- bdata[[i]]$RY[-isin_index]
     bdata[[i]]$IBOX <- bdata[[i]]$IBOX[-isin_index]
     bdata[[i]]$IBXA <- bdata[[i]]$IBXA[-isin_index]
     bdata[[i]]$IBXB <- bdata[[i]]$IBXB[-isin_index]
     bdata[[i]]$CASHFLOWS$ISIN <- bdata[[i]]$CASHFLOWS$ISIN[-cf_isin_index]
     bdata[[i]]$CASHFLOWS$CF <- bdata[[i]]$CASHFLOWS$CF[-cf_isin_index]
     bdata[[i]]$CASHFLOWS$DATE <- bdata[[i]]$CASHFLOWS$DATE[-cf_isin_index]
  }
    bdata
}

###################################################################
#               Diagnostic data checks                            #
###################################################################

# TODO: - add functionality of testscript.R

dyn_diag <- function(dynbonddata,group,fig=TRUE) {
  x <- matrix(0,ncol=length(dynbonddata),nrow=length(dynbonddata[[1]]$ISIN))
  ym <- list()
  for (i in seq(length(dynbonddata))){
    bonddata <- list()
    bonddata[[group]] <- dynbonddata[[i]]
    ym[[i]] <- diag_bond(bonddata)
    x[,i] <- diag_bond(bonddata)$y[[1]][,2]

  }
   
  if(fig)  {
    plot(x[1,],ylim=c(min(x),max(x)),type="n")
    for(i in seq(nrow(x)))  lines(x[i,])
  }
  res <- list(x=x,ym=ym)
  res
}
