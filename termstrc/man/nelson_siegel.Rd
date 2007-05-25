\name{nelson_siegel}
\alias{nelson_siegel}

\title{ Spot rate function according to Svensson }
\description{
    This function calculates the spotrates for certain maturity dates and a parameter vector
  according to Nelson and Siegel.
}
\usage{
nelson_siegel(beta, m)
}

\arguments{
  \item{beta}{ a vector of parameters \eqn{(\beta_0, \beta_1, \beta_2, \tau_1)} }
  \item{m}{ one maturity (or a vector of maturities) }
  }
}
%\details{}
 
\value{returns a vector consisting of the calculated spot rates}

\references{ Charles R. Nelson and Andrew F. Siegel (1987): Parsimonious modeling of yield curves. \empd{The Journal 
of Business}, \bold{60(4)}:473--489.}

\author{ Robert Ferstl, Josef Hayden }

\note{ For details please refer to the vignette of the package or the article of Nelson and Siegel.
}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

nelson_siegel(rep(0.01,4),1:30)
}
\keyword{  }

