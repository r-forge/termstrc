\name{aabse}
\alias{aabse}

\title{Average Absolute Mean Error}

\description{Average absolute mean error}
\usage{
aabse(actual, estimated)
}

\arguments{
  \item{actual}{vector, consisting of the observered values, e.g. bond prices}
  \item{estimated}{vector, consisting of the estimated values, e.g. bond prices}
}
\details{
  Calculation of the AABSE according to the formula:
\deqn{AABSE=\frac{1}{n}\sum{|e-\overline{e}|}}
}
%\value{}
\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canda.\emph{Technical Report No 84 Bank of Canada}}

\author{Robert Ferstl, Josef Hayden}
\seealso{\code{\link{rmse}}}

\keyword{}
