\name{estimateyieldcurve}
\alias{estimateyieldcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate Zero-coupon Yield Curves
}
\description{
Estimate Zero-coupon Yield curves assuming a certain spot rate function
}
\usage{
estimateyieldcurve(y, m, beta, objfct, grad_objfct, constraints, constrOptimOptions)
}

\arguments{
  \item{y}{yields
}
  \item{m}{maturities
}
  \item{beta}{parameter vector
}
  \item{objfct}{objective function}
  \item{grad_objfct}{grad_objfct
}
  \item{constraints}{parameter constraints
}
  \item{constrOptimOptions}{solver options
}
}
\details{internal helper function}
