\name{duration}
\alias{duration}
\title{Duration, modified Duration and Duration based Weights}
\description{
  Duration, modified duration and duration based weights}
\usage{
duration(cf_p, m_p, y)
}

\arguments{
  \item{cf_p}{cashflow matrix (including prices)}
  \item{m_p}{maturity matrix (including prices)}
  \item{y}{yield of the bond}
}
%\details{}
\value{
  The function returns a matrix with three columns, i.e. duration,
  modified duration and duration based weights. 
}
\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canda.\emph{Technical Report No 84 Bank of Canada}}

\author{Robert Ferstl, Josef Hayden }

%\seealso{  }
%\examples{}
\keyword{  }
