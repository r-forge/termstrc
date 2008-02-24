\name{splines_estim}
\alias{splines_estim}

\title{Discount Curve Estimation with McCulloch Cubic Splines}
\description{
  Discount curve estimation with the cubic splines approach by McCulloch (1975). }

\usage{
splines_estim(group,
              bonddata,
              matrange = "all")
}

\arguments{
  \item{group}{vector defining the group of bonds used for the estimation,\newline
               e.g. \code{c("GERMANY","AUSTRIA")}.}
  \item{bonddata}{a data set of bonds in list format.}
  \item{matrange}{use \code{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \code{c(lower,upper)}. }
}
\details{
 \item{\code{group}}{The first element of the vector will be used as the
  			    	  	reference country for the spread curve calculation.
					 \code{group} can be either a vector of bond groups or a scalar.}
\item{\code{bonddata}}{The package is designed to work with a certain list data structure. For more information use the function \code{str()} to explore the structure of the example data sets.}
}

\value{
  The function \code{splines_estim} returns an object of the class \code{"cubicsplines"}. The object
contains the following items (mainly lists):

  \item{group}{group of bonds (e.g. countries) used for the estimation.}
  \item{matrange}{\code{"none"} or a vector with the maturity range.}
  \item{n_group}{length of object \code{group}, i.e. the number of countries.}
  \item{knotpoints}{selected knot points for the cubic splines estimation.}
 \item{spot}{zero-coupon yield curves as object of the class \code{"spot_curves"}.}
  \item{spread}{spread curves as object of the class \code{"s_curves"}.}
  \item{forward}{forward curves as object of the class \code{"fwr_curves"}.} 
 \item{discount}{discount curves as object of the class \code{"df_curves"}.} 
 \item{cf}{cashflow matrices.}
  \item{m}{maturity matrices.}
  \item{p}{dirty prices.}
  \item{phat}{estimated bond prices.}
  \item{perrors}{pricing errors and maturities as object of the class \code{"error"}.}
  \item{y}{bond yields.}
  \item{yhat}{one list for each group with the theoretical bond yields calculated with the estimated bond prices \code{phat}.}
  \item{yerrors}{yield errors and maturities as object of the class \code{"error"}.}
  \item{alpha}{OLS coefficients of cubic splines estimation.}
  \item{regout}{OLS estimation results as object of the class \code{"lm"}.}
}
\references{
J.Huston McCulloch (1971):
Measuring the Term Structure of Interest Rates.
\emph{The Journal of Business, \bold{44} 19--31.} 

J. Huston McCulloch (1975):
The Tax-Adjusted Yield Curve.
\emph{The Journal of Finance, \bold{30} 811--830.} 
}

%\author{Robert Ferstl, Josef Hayden}

\note{ For objects of the class \code{"spot_curves"}, 
\code{"s_curves"}, \code{"df_curves"}, \code{"fwr_curves"}, \code{"error"} appropriate plot methods are offered. For objects of the list item \code{regout} standard \code{lm} methods apply. For objects of the class \code{"cubicsplines"} print, summary and plot methods are available. Another term structure estimation method is provided by the function \code{\link{nelson_estim}}.} 

\seealso{\code{\link{print.cubicsplines}}, \code{\link{summary.cubicsplines}},  \code{\link{plot.cubicsplines}},  \code{\link{nelson_estim}},  \code{\link{plot.spot_curves}},  \code{\link{plot.s_curves}},  \code{\link{plot.df_curves}},  \code{\link{plot.fwr_curves}},  \code{\link{plot.error}},  \code{\link{summary.lm}},  \code{\link{print.lm}},  \code{\link{plot.lm}}.}

\examples{
# load data set
data(eurobonds)

# define countries, for which the estimation 
# of the zero-coupon yield curves will be carried out
group <- c("GERMANY", "AUSTRIA", "ITALY")

# define data set
bonddata <- eurobonds

# set maturtiy range
matrange <- c(0, 19)  

# perform estimation
x <- splines_estim(group, bonddata, matrange)

# print the obtained parameters of the estimation
print(x)

# goodness of fit measures
summary(x)

# plot the zero-coupon yield curve for each country
plot(x,errors="none")

# plot all zero-coupon yield curves together
plot(x,multiple=TRUE,errors="none")

# spread curve splot
plot(x,ctype="spread",errors="none")

# price error plot for all countries
plot(x,ctype="none")
}

\keyword{}

