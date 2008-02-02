\name{nelson_siegel}
\alias{nelson_siegel}

\title{Spot Rate Function according to Nelson/Siegel}
\description{
    This function calculates the spot rates for certain maturity dates and a parameter vector
  according to Nelson/Siegel (1987).
}
\usage{
nelson_siegel(beta, m)
}

\arguments{
  \item{beta}{ a vector of parameters \eqn{(\beta_0, \beta_1, \beta_2, \tau_1)}.}
  \item{m}{ one maturity (or a vector of maturities).}
  }
}
\details{The spot rates according to Nelson/Siegel are defined as:
\deqn{
  s(m,\bm{b}) = \beta_0 + \beta_1\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} + \beta_2\left(\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} - \exp(-\frac{m}{\tau_1})\right).}}
 
\value{Returns a vector consisting of the calculated spot rates.}

\references{ Charles R. Nelson and Andrew F. Siegel (1987): Parsimonious modeling of yield curves. \emph{The Journal 
of Business}, \bold{60(4)}:473--489.}

%\author{ Robert Ferstl, Josef Hayden }


\examples{

nelson_siegel(rep(0.01,4),1:30)
}
\keyword{  }

