\name{duration}
\alias{duration}
\title{ Function for the calculation of the duration.}
\description{
Function for the calculation of the duration, the modified duration and the weights used for
the optimisation}
\usage{
duration(cf, m, y)
}

\arguments{
  \item{cf}{cashflow matrix }
  \item{m}{maturity matrix }
  \item{y}{yield of the bond}
}
%\details{}
\value{
  The function returns a matrix, with fhe following characteristics: 
  
  \item{column 1 }{Consists of the calculated duration}
  \item{column 2}{Consists of the calculated modified duration}
  \item{column 3}{Consists of the calculated duration based weights}
}
\references{Yield Curve Modelling at the Bank of Canda, David Bolder and David Streliski, Technical Report No 84 Bank of Canada 1999, p.42 }

\author{Robert Ferstl, Josef Hayden }

\note{For details please refer to the vignette of the package or to the reference} 


%\seealso{  }
%\examples{}
\keyword{  }
