\name{estimatezcyieldcurve}
\Rdversion{1.1}
\alias{estimatezcyieldcurve}

\title{Estimate Zero-coupon Yield Curves
}
\description{
Estimate Zero-coupon Yield curves assuming a certain spot rate function
}
\usage{
estimatezcyieldcurve(method, startparam, obj_fct, constraints, constrOptimOptions)
}

\arguments{
  \item{method}{form of the spot rate function
}
  \item{startparam}{start parameter vector 
}
  \item{obj_fct}{objective function, e.g., sum of the weighted squared price errors}
  \item{constraints}{constraints for the solver}
  \item{constrOptimOptions}{solver options
}
}
\details{Used as internal helper function 
}
