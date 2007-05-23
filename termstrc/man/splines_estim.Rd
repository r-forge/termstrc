\name{splines_estim}
\alias{splines_estim}

\title{ ~~function to do ... ~~ }
\description{
  Cubic splines estimation of zero coupon yield curve
}
\usage{
splines_estim(group, bonddata, maturity_spectrum = "all")
}

\arguments{
  \item{group}{a group of bonds from a data set, e.g. countries, rating classes}
  \item{bonddata}{a data set of bonds in list format}
  \item{maturity_spectrum}{maturity range of bonds from the data set used in estimation}
}
\details{
  bla bla
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{
Yield Curve Modelling at the Bank of Canda, David Bolder and David Streliski, Technical Report No 84 Bank of Canada 1999.
}

\author{Robert Ferstl, Josef Hayden}

\seealso{for another estimation method see \code{\link{nelson_estim}}}

\examples{
group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- "all"

x <- splines_estim(group, bonddata, matrange)
print(x)
summary(x)
plot(x)
}

\keyword{}

