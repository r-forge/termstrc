\name{aabse}
\alias{aabse}

\title{Average Absolute Mean Error}

\description{Calculation of the average absolute mean error (AABSE).}
\usage{
aabse(actual, estimated)
}

\arguments{
  \item{actual}{vector, consisting of the observed values.}
  \item{estimated}{vector, consisting of the estimated values.}
}


\details{
  Calculation of the AABSE according to the formula:
\deqn{AABSE=frac{1}{m}\bm{\iota}'|\bm{\epsilon}|.}
Whereas \eqn{\epsilon} is the vector of the yield or price errors of the bonds. 

}
%\value{}
%\references{}

%\author{Robert Ferstl, Josef Hayden}
\seealso{\code{\link{rmse}}}

\keyword{}
