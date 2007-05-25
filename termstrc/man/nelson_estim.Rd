\name{nelson_estim}
\alias{nelson_estim}

\title{Function for term structure and credit spread estimation }
\description{
  Main function for term structure and credit spread estimation based on the approaches of Nelson and Siegel and Svensson.
}
\usage{
nelson_estim(group, bonddata, matrange = "all", method = "Nelson/Siegel", fit = "prices", weights = "none", startparam, control = list(eval.max = 1000))
}


\arguments{
  \item{group}{vector defining the group of bonds used for the estimation, e.g. \code{c("GERMANY","AUSTRIA") }}
  \item{bonddata}{list with bonddata}
  \item{maturity_spectrum}{ use \code{"all"} for no restrictions, or restrict the bonds used for the estimation \code{c(lower,uppper)} }
  \item{method}{\code{"Nelson/Siegel"} or \code{"Svensson"} }
  \item{fit}{Use \code{"prices"} (\code{"yields"}) for minimising the squared price (yield) deviation}
  \item{weights}{If a weighted minimisation of the squared price deviation is required, use \code{"duration"}, otherwise \code{"none"} }
  \item{startparam}{ matrix of startparameters, for the Nelson/Siegel (Svensson) method 4 (6) parameters for each each group are required (one row for one group).}
  \item{control}{list of control parameters for the function \code{nlminb}}
}



\details{
 
 
 \item{\code{group}}{The first element of the vector will be used as the reference country for the credit spread estimation.
The group can be either a vector of groups or a scalar.}
 
\item{\code{bonddata}}{NOTE: The package is tailored to the included data sets. Therefore the structure and the naming convention of other used data sets has to be identical.  
 Use the function \code{str()} to explore the structure of the provided datasets.}
 
 \item{\code{startparam}}{For the estimation a matrix of start parameters must be provided. The estimation is very sensitive to the start parameters, due to multiple
 local minima of the objective functions. Therefore a great effort should be invested in finding a good start parameterset. }
 
 \item{\code{control}}{Please refer to the documentation of the function \code{mlimnb}.}

 
 }
\value{

The function \code{nelson_estim} returns a list with the following elements or sub-lists:

  \item{group}
  \item{matspectrum }{Includes the chosen maturity range}
  \item{method}{Includes the chosen estimation method}
  \item{fit}{Includes the chosen objective function}
  \item{weights}{Includes the type of optimisation, e.g "none" for a non-weighted optimisation}
  \item{n_group}{the number of countries (rating categories) used for the optimisation}
  \item{zcy_curves}{Includes the values for plotting the estimated zero-coupon yield curves}
  \item{scuves} {Includes the values for plotting the spread curves}
  \item{cf}{the cashflow matrix for all specified groups}
  \item{m}{the maturity matrix for all specified groups}
  \item{duration}{includes the duration, weighted duration and duration based weights}
  \item{p}{the dirty prices for all specified groups}
  \item{phat}{the estimated prices for all specified groups}
  \item{y}{the calculated yields for all bonds of the specified countries}
  \item{yhat}{the theoretical yields calculated on the basis of \code{phat}}
  \item{opt_result}{the optimal parameter vector for the specified groups according to the chosen estimation approach}
   
}



\references{ 
Charles R. Nelson and Andrew F. Siegel (1987): 
Parsimonious modeling of yield curves.
\empd{The Journal of Business}, \bold{60(4)}:473--489.

Lars E.O. Svensson (1994): 
Estimating and interpreting forward interest rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}. 
 
David Bolder and David Streliski (1999):
Yield Curve Modelling at the Bank of Canada.
\emph{Technical Report No 84 Bank of Canada} 

 }

\author{Robert Ferstl, Josef Hayden}


%\note{ }


%\seealso{\code{\link{help}}, ~~~ }
\examples{ 
demo(euro01)
}


\keyword{  }

