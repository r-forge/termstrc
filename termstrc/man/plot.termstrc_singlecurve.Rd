\name{plot.termstrc_singlecurve}
\alias{plot.termstrc_singlecurve}

\title{ S3 plot method}
\description{
 S3 plot method
}
\usage{
plot.termstrc_singlecurve(x, maturity_min = 1, maturity_max = 15, spread_curves = TRUE,...)
}
\arguments{
  \item{x}{ object of the class \code{termstrc_singlecurve}}
  \item{maturity_min}{ defines the minimum maturity used for the plots }
  \item{maturity_max}{ defines the minimum maturity used for the plots }
  \item{spread_curves}{ If \code{TRUE} spreadcurves are plotted }
  \item{...}{optional arguments}	

}

\details{A plot of the estimated zero coupon yield curve for each specified group of bonds (e.g 
countries, ratingclasses) is offered. In Addition the spread curves are plotted, if desired. }
 
%\value{ }
%\references{ ~put references to the literature/web site here ~ }
\author{Robert Ferstl, Josef Hayden}

%\note{ }

\seealso{ \code{\link{print.termstrc_singlecurve}},\code{\link{summary.termstrc_singlecurve}}}
%\examples{}
\keyword{  }

