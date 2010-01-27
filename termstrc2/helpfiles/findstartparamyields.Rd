\name{findstartparamyields}
\alias{findstartparamyields}

\title{Find Globally Optimal Start Parameters}
\description{Start parameter search routine for term structure estimation based on a yield data set. 
The algorithm searches for the parameters over a grid spanned over tau1 (tau2). 
}


\usage{
findstartparamyields(y, m, method, deltatau = 0.1)
}
\arguments{
  \item{y}{yields
}
  \item{m}{maturities
}
  \item{method}{type of spot rate function
}
  \item{deltatau}{stepsize of the grid 
}
}
\details{Used as internal helper function
}
\value{Returns an object of the class \command{"spsearch"}, which
  includes the startparameters and details concerning the optimization. 
}