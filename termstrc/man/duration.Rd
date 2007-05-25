\name{duration}
\alias{duration}
\title{Calculation of the duration, modified duration and duration based weights.}
\description{
Function for the calculation of the duration, the modified duration and duration based weights}
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
\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canda.\emph{Technical Report No 84 Bank of Canada}}

\author{Robert Ferstl, Josef Hayden }

\note{For details please refer to the vignette of the package or to the reference} 


%\seealso{  }
%\examples{}
\keyword{  }
