\name{bond_yields}
\alias{bond_yields}

\title{Bond yield calculation}

\description{Function for the calculation of bond yields.}

\usage{
bond_yields(cashflows, m, tol = 1e-10)
}

\arguments{
  \item{cashflows}{ matrix with the bonds cashflows }
  \item{m}{maturity matrix}
  \item{tol}{desired accuracy for function \code{uniroot} }
}
%\details{}

\value{
matrix with the bond yields and the associated maturity 
}

\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canda.\emph{Technical Report No 84 Bank of Canada}}

\author{Robert Ferstl, Josef Hayden }
%\note{}

\seealso{ \code{\link{uniroot}}}
\examples{

}
\keyword{ }
