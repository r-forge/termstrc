\name{bond_prices}
\alias{bond_prices}

\title{Bond price calculation }
\description{
  Function for the calculation of bond prices according to the chosen approach 
  (Nelson and Siegel or Svensson) based on the cashflows and maturities of the bonds. 
  
}
\usage{
bond_prices(method = "Nelson/Siegel", beta, m, cf)
}

\arguments{
  \item{method}{ defines the desired method,  \code{"Nelson/Siegel"} for the Nelson/Siegel approach
  		or \code{"Svensson"} for the Svensson approach}
  \item{beta}{parameter vector, is linked to the chosen approach}
  \item{m}{maturities matrix, consists
  of the maturity dates which are appendant to the cashflows of the bonds   }
  \item{cf}{ cashfows matrix}
 }
%\details{}
\value{
  
  The function bond prices returns a list with:
  
  \item{spot_rates}{spot rates}
  \item{discount_factors}{discount factors}
  \item{bond_prices}{bond prices}
}

\references{\emph{Yield Curve Modelling at the Bank of Canda, David Bolder and David Streliski}, Technical Report No 84 Bank of Canada 1999, p.39 }

\author{ Robert Ferstl, Josef Hayden }
\note{  


}
\seealso{ \code{\link{svensson}}, \code{\link{nelson_siegel}}}

%\examples{}

\keyword{  }