\name{loss_function}
\alias{loss_function}

\title{Loss function used for the term structure estimation }
\description{
  The loss function defines the objective function used for the optimisation. 
  
  
}
\usage{
loss_function(p, phat, omega, weights)
}

\arguments{
  \item{p}{ The observed vector of prices, yields }
  \item{phat}{ The estimated vector of prices, yields }
  \item{omega}{ Weights }
  \item{weights}{ if \code{"none"} the squared deviation of \code{p} and \code{phat} is not weighted, use \code{"duration"} 
  					for a duration weighted optimisation}
}
\details{
 For details please refer to the vignette of the package or to reference.
}
%\value{}
\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canda.\emph{Technical Report No 84 Bank of Canada}}
\author{Robert Ferstl, Josef Hayden}
%\note{ ~~further notes~~}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{}
\keyword{}
