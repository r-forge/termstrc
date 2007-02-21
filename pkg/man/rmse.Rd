\name{rmse}
\alias{rmse}
\title{Root mean-square error}
\description{Root mean square error}

\usage{
rmse(actual, estimated)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{actual}{vector, consisting of the observered values, e.g. observed bond prices (yields)}
  \item{estimated}{vector, consisting of the estimated values, e.g. estimated bond prices (yields) }
}
\details{
  Calculation of the RMSE according to the formula:
\deqn{RMSE=\sqrt{\frac{1}{n}\sum{(e-\overline{e})^2}}}

}
\references{Yield Curve Modelling at the Bank of Canda, David Bolder and David Streliski, Technical Report No 84 Bank of Canada 1999, p.25 }

\author{Robert Ferstl, Josef Hayden}

\seealso{\code{\link{aabse}}}

\keyword{}
