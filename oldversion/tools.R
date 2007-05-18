###################################################################
#            cashflows matrix for a given country                 #
###################################################################

create_cashflows_matrix <- function(country,include_price=FALSE) {
  
  n_of_cf <- summary(as.factor(country$CASHFLOWS$ISIN))
  n_of_bonds <- length(n_of_cf)
  max_cf <- max(n_of_cf)
  pos_cf <- c(0,cumsum(n_of_cf))
  
  #calculate cashflow matrix:
  #the number of rows of the matrix = number of 
  #cashflows of the bond with the maximum of cashflows
  #all missing elements of the matrix are filled up with zeros 
  	
  CASHFLOWMATRIX <-
    mapply(function(i) c(country$CASHFLOWS$CF[(pos_cf[i]+1):pos_cf[i+1]],
                         rep(0,max_cf-n_of_cf[i])),
           1:n_of_bonds)
  
  if (include_price == TRUE) {CASHFLOWMATRIX <- rbind(-(country[["PRICE"]]
      +country[["ACCRUED"]]),CASHFLOWMATRIX)}              
  colnames(CASHFLOWMATRIX) <- country$ISIN
  CASHFLOWMATRIX
}

###################################################################
#              maturities matrix for a given country              #
###################################################################

create_maturities_matrix <-
  function(country,include_price=FALSE) {

  n_of_cf <- summary(as.factor(country$CASHFLOWS$ISIN))
  n_of_bonds <- length(n_of_cf)
  max_cf <- max(n_of_cf)
  pos_cf <- c(0,cumsum(n_of_cf))
  year_diff <- as.numeric(difftime(as.Date(country$CASHFLOW$DATE),
  				as.Date(country$TODAY),units="days"))/365
 
  #calculate maturity matrix
  #calculate cashflow matrix:
  #the number of rows of the matrix = number of 
  #maturity dates of the bond with the longest maturity
  #all missing elements of the matrix are filled up with zeros 
  
  MATURITYMATRIX <-
     mapply(function(i) c(year_diff[(pos_cf[i]+1):pos_cf[i+1]],
                     rep(0,max_cf-n_of_cf[i])),
            1:n_of_bonds)
  
  if (include_price == TRUE) {MATURITYMATRIX <- rbind(rep(0,n_of_bonds),
                           MATURITYMATRIX)}  
  colnames(MATURITYMATRIX) <- country$ISIN
  MATURITYMATRIX
}

###################################################################
#              bond yields for a given country                    #
###################################################################

bond_yields <- function(cashflows, m,tol=1e-10) {

  # convert input data to matrices, if necessary
  if (!is.matrix(cashflows))
    cashflows <- as.matrix(cashflows)
  if (!is.matrix(m))
    m < -as.matrix(m)

  # create empty bond yields matrix in appropriate size
  bondyields<-matrix(0, nrow=ncol(cashflows), ncol=2)                                                                                                                     

  # put maximum of m of every bond into first column of bond yields matrix
  bondyields[,1] <- apply(m, 2, max)

  # traverse list of bonds
  for (i in seq(ncol(cashflows))) {

  # calculate bond price with yield 
  yield_function<-function(y) {
       t(cashflows[,i])%*%exp(-m[,i]*y)
    }

  # calculate roots and put result into second column of bond yields matrix
    
  bondyields[i,2]<-uniroot(yield_function, c(0, 1), tol = tol,maxiter=3000)$root
    
  }

  # return calculated bond yields matrix
  rownames(bondyields) <- colnames(cashflows)
  colnames(bondyields) <- c("Maturity","Yield")
  bondyields
}

###################################################################
#                      maturity range                             #
###################################################################

maturity_range <- 
  function(bonddata,lower,upper) {

m <- lapply(bonddata,create_maturities_matrix)
colmax <- function(m) apply(m,2,max)
m_max <- lapply(m,colmax)

bonds_in_range <- function(country) names(country[which(
                country>lower & country<upper)])

# list with ISINs of bonds in range

isins_range <- lapply(m_max,bonds_in_range)

index_set <- which(unlist(lapply(isins_range,length))>0)

# list with positions of bonds in isins_range

isins_range_pos <- list()

for (i in 1:length(bonddata)) {
        isins_range_pos[[i]] <- which(bonddata[[i]]$ISIN %in% 
            isins_range[[i]])
    }
names(isins_range_pos) <- names(bonddata)

# first part of bonddata for filtering

if(length(bonddata[[1]])>8) N=8 else N=6

print(N)

first <- function(lst) lst[c(1:N)]
filtered <- lapply(bonddata,first)

bonddata_range <- list()
for (i in 1:length(bonddata)) {
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
for (i in 1:length(bonddata)) {
isins_range_pos[[i]] <- which(bonddata[[i]][["CASHFLOWS"]]
                        [["ISIN"]]%in%isins_range[[i]])
}
names(isins_range_pos) <- names(bonddata)

for (i in 1:length(bonddata)) {
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

for (i in 1:length(bonddata)) {
bonddata_range[[i]][["TODAY"]] <- bonddata[[i]][["TODAY"]]
}

# delete countries where no bonds are available

bonddata_range <- bonddata_range[index_set]
bonddata_range
}

###################################################################
#                      Root mean squared error                    #
###################################################################

rmse <-
function (actual,estimated) {
			e <- actual - estimated
			rmse <- sqrt(1/length(e)*sum((e-mean(e))^2))
      rmse				
      }

###################################################################
#                 Average absolute error                          #
###################################################################
      
aabse <-
function (actual,estimated){
			e <- actual - estimated	
      aabse <- 1/length(e)*sum(abs(e-mean(e)))
		}      

###################################################################
#                 Duration                                        #
###################################################################

# calculates duration, modified duration and weights for optimization

# cf ... cashflows (dirty price included)
# m ... maturities (zero included)
# y ... bond yields
duration <-
function (cf,m,y) {
       y <- matrix(rep(y,nrow(m)),ncol=ncol(m),byrow=TRUE)
       d <- apply(cf*m*exp(-y*m),2,sum)/-cf[1,]
       md <- d/(1+y[1,])
       omega <- (1/md)*sum(1/md)
       cbind(d,md,omega)
    }

# example

#cf = matrix(c(-103, 5, 5, 105,-102,3,3,103),ncol=2)
#m = matrix(c(0,1,2,3,0,0.5,1.5,2),ncol=2)
#y = bond_yields(cf,m)
#y = matrix(y[,2],ncol=2)

#duration(cf,m,y)

###################################################################
#                 Spotrate calculation                            #
###################################################################
#calculate spotrate according to chosen approach, optimal parametervector
#and maturity-vector.
 
 srates <- function(method,beta,m){
 		 
  func<- switch(method,
 				"Nelson/Siegel" = nelson_siegel(beta,m),
 				"Svensson" = svensson(beta,m))
 func
 	}

  				
  
 
 							








