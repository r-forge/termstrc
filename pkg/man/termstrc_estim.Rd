\name{termstrc_estim}
\alias{termstrc_estim}

\title{Main function for term structure and credit spread estimation }
\description{
  Main function for term structure and credit spread estimation
}
\usage{
termstrc_estim(countries, bonddata, maturity_spectrum = "all", method = "Nelson/Siegel", fit = "prices", weights = "none", startparam, control = list(eval.max = 1000))
}


\arguments{
  \item{countries}{vector defining the group of bonds used for the estimation, e.g. \code{c("GERMANY","AUSTRIA") }}
  \item{bonddata}{list with bonddata}
  \item{maturity_spectrum}{ use \code{"all"} for no restrictions, or restrict the bonds used for the estimation \code{c(lower,uppper)} }
  \item{method}{\code{"Nelson/Siegel"} or \code{"Svensson"} }
  \item{fit}{Use \code{"prices"} (\code{"yields"}) for minimising the squared price (yield) deviation}
  \item{weights}{If a weighted minimisation of the squared price(yield) deviation is required, use \code{"duration"}, otherwise \code{"none"} }
  \item{startparam}{ matrix of startparameters, for the Nelson/Siegel (Svensson) method 4 (6) parameters for each each country/rating category are required (one row for one country).}
  \item{control}{list of control parameters for the function \code{nlminb}}
}



\details{
 
 
 \item{\code{countries}}{The first element of the vector will be used as the reference country for the credit spread estimation.
 The estimation is designed for the estimation of more countries, rating categories. Is only one country specified an error message will occur.}
 
\item{\code{bonddata}}{NOTE: The package is tailored to the included data sets. Therefore the structure and the naming convention of other used data sets has to be identical.  
 Use the function \code{str()} to explore the structure of the provided datasets.}
 
 \item{\code{startparam}}{For the estimation a matrix of start parameters must be provided. The estimation is very sensitive to the start parameters, due to multiple
 local minima of the objective functions. Therefore a great effort should be invested in finding a good start parameterset. }
 
 \item{\code{control}}{Please refer to the documentation of the function \code{mlimnb}.}

 
 }
\value{

The function \code{termstrc} returns a list with the following elements or sub-lists:

  \item{maturity_spectrum }{Includes the chosen maturity spectrum}
  \item{method}{Includes the chosen estimation method}
  \item{fit}{Includes the chosen objective function}
  \item{weights}{Includes the type of optimisation, e.g "none" for a non-weighted optimisation}
  \item{n_countries}{the number of countries (rating categories) used for the optimisation}
  \item{cashflows}{the cashflow matrix for all specified countries}
  \item{maturities}{the maturity matrix for all specified countries}
  \item{dirty_prices}{the dirty prices for all specified countries}
  \item{estimated_prices}{the estimated prices for all specified countries}
  \item{yields}{the calculated yields for all bonds of the specified countries}
  \item{opt_result}{the optimal parameter vector for the specified countries according to the chosen estimation approach}
  \item{spotrates}{the spotrates calculated with the optimal parameter vector}
  
}



\references{  Charles R. Nelson and Andrew F. Siegel. Parsimonious modeling of yield curves. The Journal of Business, 60(4):473Ð489, Oct 1987.

Lars E.O. Svensson. Estimating and interpreting forward interest rates: Sweden 1992 -1994. Technical Reports 4871, National Bureau of Economic Research, Inc, Sep 1994. 

Yield Curve Modelling at the Bank of Canda, David Bolder and David Streliski, Technical Report No 84 Bank of Canada 1999.
 }

\author{Robert Ferstl, Josef Hayden}


%\note{ }


%\seealso{\code{\link{help}}, ~~~ }
\examples{ 
demo(euro01)
}


\keyword{  }

