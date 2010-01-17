\name{estim_nss.zeroyields}
\alias{estim_nss.zeroyields}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estim_nss.zeroyields(obj, method = "ns", deltatau = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{deltatau}{
%%     ~~Describe \code{deltatau} here~~
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
function (obj, method = "ns", deltatau = 1)
  {
    if(method=="ns"){
          spsearch <- findstartparamyields(obj$yields[1,],obj$maturities, method, deltatau)
          startparam <- spsearch$startparam
          objfct <- objfct_ns
          grad_objfct <- grad_objfct_ns
          ui <- rbind(c(1,0,0,0),  # beta0 > 0
                      c(1,1,0,0),  # beta0 + beta1 > 0
                      c(0,0,0,1),  # tau1 > 0
                      c(0,0,0,-1)) # tau1 < 30
          ci <- c(0,0,0,-30)
    }

    if(method=="sv"){
      
          spsearch <- findstartparamyields(obj$yields[1,],obj$maturities, method, deltatau)
          startparam <- spsearch$startparam
          
          objfct <- objfct_sv
          grad_objfct <- grad_objfct_sv
          ui <- rbind(c(1,0,0,0,0,0),  # beta0 > 0
                      c(1,1,0,0,0,0),  # beta0 + beta1 > 0
                      c(0,0,0,1,0,0),  # tau1 > 0
                      c(0,0,0,-1,0,0), # tau1 < 30
                      c(0,0,0,0,0,1),  # tau2 > 0
                      c(0,0,0,0,0,-1)) # tau2 < 30
          ci <- c(0,0,0,-30,0,-30)
    }

    optresult <- list()
    optparam <- matrix(nrow = nrow(obj$yields), ncol = length(startparam)) 
    for (i in 1:nrow(obj$yields)){
    
      yields <- obj$yields[i,]

      if (i<2)
        beta <- startparam
      else
        beta <- optparam[i-1,]
        
      optresult[[i]] <- estimateyieldcurve(yields, obj$maturities, beta, objfct, grad_objfct, ui, ci)
      optparam[i,] <- optresult[[i]]$par
    }
    colnames(optparam) <- switch(method,
          "ns"=c("beta_0","beta_1","beta_2","tau_1"),
          "sv"=c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")) 
    result <- list(optparam = optparam, optresult = optresult, method = method,
                   maturities = obj$maturities, dates = obj$dates, spsearch = spsearch)
    class(result) <- "termstrc_yields"
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
