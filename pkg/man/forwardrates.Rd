\name{forwardrates}
\alias{forwardrates}

\title{Forward rate calculation }
\description{Calculate forward rates according to the Nelson/Siegel or Svensson approach.}
\usage{
forwardrates(method, beta, m)
}

\arguments{
  \item{method}{Specify the method used for the forward rate calculation (/code{Nelson/Siegel} or /code{Svensson} }
  \item{beta}{Parameters corresponding to the chosen method }
  \item{m}{maturity vector}
}
\details{}
\value{The function returns a vector with the calculated forward rates.}
  
%\references{  }
\author{ Robert Ferstl, Josef Hayden }
%\note{ }
\seealso{\code{\link{fwr_ns}}, \code{\link{fwr_sv}}}
%\examples{}

\keyword{  }
