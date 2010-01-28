\name{estim_nss.zeroyields}
\alias{estim_nss.zeroyields}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{S3 estim_nss Method
}
\description{The method performs an iterative term structure estimation procedure
  on a dynamic yield data set of the class
  \command{"zeroyields"}. Available methods are Nelson/Siegel,
  Diebold/Li and (adjusted) Svensson. 
}
\usage{
\method{estim_nss}{zeroyields}(dataset, method = "ns", lambda = 0.0609 * 12, deltatau = 1, optimtype = "allglobal", constrOptimOptions = list(control = list(), outer.iterations = 200, outer.eps = 1e-04), ...)
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
  \item{deltatau}{stepsize for parameter grid used for globally
    optimal parameter search
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

\value{ The method returns an object of the class
  \command{"dyntermstrc_yields"}.There are print, plot and summary
  method available. 
}

\examples{## Run: demos(zeroyields) }