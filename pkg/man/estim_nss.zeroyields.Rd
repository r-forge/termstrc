\name{estim_nss.zeroyields}
\alias{estim_nss.zeroyields}

\title{S3 estim_nss Method
}
\description{The method performs an iterative term structure estimation procedure
  on a dynamic yield data set of the class
  \command{"zeroyields"}. Available methods are Nelson/Siegel,
  Diebold/Li and (adjusted) Svensson. 
}
\usage{
\method{estim_nss}{zeroyields}(dataset, method = "ns", lambda = 0.0609 *12, tauconstr = NULL, optimtype = "allglobal", constrOptimOptions = list(control = list(), outer.iterations = 200, outer.eps = 1e-04), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{dynamic bond data set of the class \command{"zeroyields"}
}
  \item{method}{\command{"ns"} for Nelson/Siegel (default), \command{"dl"} for
    Diebold/Li, \command{"sv"} for Svensson or \command{"asv"} for adjusted Svensson.
}
  \item{lambda}{parameter on a yearly time scale with fixed value for \command{"dl"} spot rate
    function (default: 0.0609*12)
}
\item{tauconstr}{This is  vector with parameters for the grid search procedure containing:

  For parametrizations except Diebold/Li, a grid search for the tau-parameter is performed. The parameters must lie within the following bounds.
  
  lower bound < [tau_1, tau_2] < upper bound
  
  The width of the grid is given by gridsize.

  tau_2 - tau_1 > taudistance
  
  (upper bound, lower bound, gridsize, tau distance)
}
  \item{optimtype}{use \command{"firstglobal"} for an inital search for
  globally optimal start parameters or \command{"allglobal"} for a
  search at every iteration. 
}
  \item{constrOptimOptions}{list with solver control parameters
    (default: control=list(), outer.interations=30,
    outer.eps.=1e-04). For further documentation please refer to
    \command{\link{optim}}
}
  \item{\dots}{further arguments
}
}

\value{The method returns an object of the class
  \command{"dyntermstrc_yields"}. There are print, plot and summary
  method available. 
}

\references{
Michiel De Pooter (2007): Examining the Nelson-Siegel Class of Term Structure Models: In-Sample Fit versus Out-of-Sample Forecasting Performance.

F.X. Diebold and C. Li: Forecasting the Term Structure of Government Bond Yields. \emph{Journal of Econometrics}, \bold{130}:337--364.

Charles R. Nelson and Andrew F. Siegel (1987): Parsimonious Modeling of Yield Curves. \emph{The Journal of Business}, \bold{60(4)}:473--489.

Lars E.O. Svensson (1994): Estimating and Interpreting Forward Interest Rates: Sweden 1992 -1994. \emph{Technical Reports 4871, National Bureau of Economic Research}.

}

\examples{
## Load yields matrix from csv
data(zyields)
x <- zyields

maturities <- c(1/12,3/12,6/12,9/12,1:12)
yields <- as.matrix(x[,2:ncol(x)])
dates <- as.Date(x[,1],format="%d.%m.%Y")

## Call class constructor
datazeroyields <- zeroyields(maturities, yields, dates)

## Estimate Nelson/Siegel model
ns_res <- estim_nss(datazeroyields, "ns", tauconstr = c(0.1, 6, 0.1))
print(ns_res)
summary(ns_res)

## Plot parameter time series
plot(param(ns_res))

## Plot yield curves in 3D
plot(ns_res)

## Plot start parameter search for t=1
plot(ns_res$spsearch[[1]])

## For more examples
## Run: demo(zeroyields)}
