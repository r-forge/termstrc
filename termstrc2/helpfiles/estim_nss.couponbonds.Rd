\name{estim_nss.couponbonds}
\alias{estim_nss.couponbonds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estim_nss.couponbonds(bonddata, group, matrange = "all", method = "ns", startparam = NULL, lambda = 0.0609 * 12, deltatau = 0.1, control = list(), outer.iterations = 30, outer.eps = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bonddata}{
%%     ~~Describe \code{bonddata} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{matrange}{
%%     ~~Describe \code{matrange} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{startparam}{
%%     ~~Describe \code{startparam} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{deltatau}{
%%     ~~Describe \code{deltatau} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{outer.iterations}{
%%     ~~Describe \code{outer.iterations} here~~
}
  \item{outer.eps}{
%%     ~~Describe \code{outer.eps} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(bonddata,                  # dataset (static)
                                  group,                     # names of countries for estimation c("Country 1", "Country 2", ...)
                                  matrange="all",            # maturity range in years c(min, max) 
                                  method="ns",
                                  startparam=NULL,           # startparameter matrix with columns c("beta0","beta1","beta2","tau1","beta3","tau2")
                                        # otherwise globally optimal parameters are searched automatically
                                  lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                                  deltatau=0.1,              # interval for parameter grid
                                  control=list(),            # options or optim() 
                                  outer.iterations = 30,     # options for constrOptim()
                                  outer.eps = 1e-04
           ) {

  ## data preprocessing
  prepro <- prepro_bond(group=group,bonddata=bonddata,matrange=matrange)

  n_group=prepro$n_group
  sgroup=prepro$sgroup
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  duration=prepro$duration
 
  ## automatically determine globally optimal start parameters
  spsearch <- list()
  length(spsearch) <- n_group
  
  if(is.null(startparam)){
    startparam <- matrix(ncol = 6, nrow = n_group)
    
    colnames(startparam) <- c("beta0","beta1","beta2","tau1","beta3","tau2")
    
    if (method == "dl") {startparam <- startparam[,1:3, drop=FALSE]}
    if (method == "ns") {startparam <- startparam[,1:4, drop=FALSE]}
    
    for (k in sgroup){
      print(paste("Searching startparameters for ", group[k]))
      spsearch[[k]] <- findstartparambonds(p[[k]],m[[k]],cf[[k]], duration[[k]][,3],
                                            method, deltatau)
      startparam[k,] <- spsearch[[k]]$startparam 
      print(startparam[k,])
    }
  }

  rownames(startparam) <- group
  
  ## objective function (weighted price error minimization) 
  obj_fct <- function(b) {
     loss_function(p[[k]],
     	bond_prices(method,b,m[[k]],cf[[k]],lambda)$bond_prices,duration[[k]][,3])}
                  
  if(method=="dl"){
    ui <- rbind(c(1,0,0),               # beta0 > 0
                c(1,1,0))               # beta0 + beta1 > 0
    ci <- c(0,0)
   }
 
  if(method=="ns"){
    ui <- rbind(c(1,0,0,0),             # beta0 > 0
                c(1,1,0,0),             # beta0 + beta1 > 0
                c(0,0,0,1),             # tau1 > 0
                c(0,0,0,-1))            # tau1 < 30
    ci <- c(0,0,0,-30)
    }

   if(method=="sv"){
     ui <- rbind(c(1,0,0,0,0,0),        # beta0 > 0
                 c(1,1,0,0,0,0),        # beta0 + beta1 > 0
                 c(0,0,0,1,0,0),        # tau1 > 0
                 c(0,0,0,-1,0,0),       # tau1 < 30
                 c(0,0,0,0,0,1),        # tau2 > 0
                 c(0,0,0,0,0,-1))       # tau2 < 30
     ci <- c(0,0,0,-30,0,-30)
    }

  ## calculate optimal parameter vectors
  opt_result <- list()

  for (k in sgroup){
    opt_result[[k]] <- constrOptim(theta = startparam[k,],
                               f = obj_fct,
                               grad = NULL,
                               ui = ui,
                               ci = ci,
                               mu = 1e-04,
                               control = control,
                               method = "Nelder-Mead",
                               outer.iterations = outer.iterations,
                               outer.eps = outer.eps)              
  }

  ## data post processing 
  postpro <- postpro_bond(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,method,lambda)
  
  ## return list of results 
  result <- list(group=group,                   # e.g. countries, rating classes
                 matrange=matrange,             # maturity range of bonds
                 method=method,                 # estimation method
                 startparam=startparam,         # calculated startparameters
                 n_group=n_group,               # number of groups,
                 spsearch = spsearch,            # detailed data from start param search
                 spot=postpro$zcy_curves,       # zero coupon yield curves
                 spread=postpro$s_curves,       # spread curves
                 forward=postpro$fwr_curves,    # forward rate curves
                 discount=postpro$df_curves,    # discount factor curves
                 expoints=postpro$expoints,     # extrapolation points
       		 cf=cf,                         # cashflow matrix
                 m=m,                           # maturity matrix
                 duration=duration,             # duration, modified duration, weights
                 p=p,                           # dirty prices        
                 phat=postpro$phat,             # estimated dirty prices         
                 perrors=postpro$perrors,       # price errors
                 ac=ac,                         # accrued interest
                 y=y,                           # maturities and yields
                 yhat=postpro$yhat,             # estimated yields
                 yerrors=postpro$yerrors,       # yield errors
                 opt_result=opt_result          # optimisation results           
                 )
              
  for ( i in 6:length(result)) names(result[[i]]) <- group
  class(result) <- "termstrc_ns"
  result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
