\name{create_cashflows_matrix}
\alias{create_cashflows_matrix}

\title{ Function for the creation of the cashflow matrix }
\description{
  Creates a matrix of cashflows for a specified country with bonds 
  in columns. The number of rows is the number of cashflows for the bond with the longest maturity.

}
\usage{
create_cashflows_matrix(country, include_price = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{country}{country name}
  \item{include_price}{If \code{TRUE} the dirty price is included }
}


\value{matrix}
%\references{ ~put references to the literature/web site here ~ }
\author{Robert Ferstl, Josef Hayden}
%\note{ }
\seealso{ \code{\link{create_maturities_matrix}} }
%\examples{}
\keyword{ }

