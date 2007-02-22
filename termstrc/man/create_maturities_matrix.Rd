\name{create_maturities_matrix}
\alias{create_maturities_matrix}

\title{ Function for the creation of the maturity matrix  }
\description{
  Creates a matrix of maturities for a specified country with bonds 
  in columns. The number of rows is the number of cashflows for the bond with the longest maturity.

}

\usage{
create_maturities_matrix(country, include_price = FALSE)
}

\arguments{
  \item{country}{country name}
  \item{include_price}{ If TRUE the dirty price is included }
}


\value{matrix
}

%\references{ ~put references to the literature/web site here ~ }
\author{Robert Ferstl, Josef Hayden}

%\note{ ~~further notes~~ 

\seealso{\code{\link{create_cashflows_matrix}}}

%\examples{}
\keyword{ }
