\name{bond_yields}
\alias{bond_yields}

\title{ Function for the calculation of bond yields}

\description{Function for the calculation of the bond yields. }

\usage{
bond_yields(cashflows, m, tol = 1e-10)
}

\arguments{
  \item{cashflows}{ matrix with the bonds cashflows }
  \item{m}{ maturity matrix}
  \item{tol}{ tol = desired accuracy for function \code{uniroot} }
}
%\details{}

\value{
matrix with the bond yields and the associated maturity 
}

\references{Yield Curve Modelling at the Bank of Canda, David Bolder and David Streliski, Technical Report No 84 Bank of Canada 1999, p.41 }

\author{ Robert Ferstl, Josef Hayden }
%\note{}

\seealso{ \code{\link{uniroot}}}
\examples{

}
\keyword{ }
