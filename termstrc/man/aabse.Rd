\name{aabse}
\alias{aabse}

\title{Average absolute mean error}

\description{Average absolute mean error}
\usage{
aabse(actual, estimated)
}

\arguments{
  \item{actual}{vector, consisting of the observered values, e.g. observed bond prices }
  \item{estimated}{vector, consisting of the estimated values, e.g. estimated bond prices}
}
\details{
  Calculation of the AABSE according to the formula:
\deqn{AABSE=\frac{1}{n}\sum{|e-\overline{e}|}}
}
%\value{}
\references{Yield Curve Modelling at the Bank of Canda, David Bolder and David Streliski, Technical Report No 84 Bank of Canada 1999, p.25 }

\author{Robert Ferstl, Josef Hayden}
\seealso{\code{\link{rmse}}}

\keyword{}
