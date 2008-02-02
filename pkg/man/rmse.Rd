\name{rmse}
\alias{rmse}
\title{Root Mean Squared Error}
\description{Calulates the root mean squared error (RMSE).}

\usage{
rmse(actual, estimated)
}

\arguments{
  \item{actual}{vector, consisting of the observed values.}
  \item{estimated}{vector, consisting of the estimated values.}
}
\details{
  Calculation of the RMSE according to the formula:
\deqn{ \mbox{RMSE}=\sqrt{\frac{1}{m}\bm{\iota}'\epsilon^2} }.
Whereas \eqn{\epsilon} is the vector of the yield or price errors of the bonds. 
}

\seealso{\code{\link{aabse}}}

\keyword{}
