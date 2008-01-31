\name{plot.s_curves}
\alias{plot.s_curves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.s_curves(x, xlim = c(range(mapply(function(i) range(x[[i]][, 1]), seq(x)))), ylim = c(range(mapply(function(i) range(x[[i]][, 2]), seq(x)))) * 10000, expoints = NULL, xlab = "Maturity (years)", ylab = "Spread (basis points)", lwd = 2, lty = 1, main = "Spread curves", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{xlim}{ ~~Describe \code{xlim} here~~ }
  \item{ylim}{ ~~Describe \code{ylim} here~~ }
  \item{expoints}{ ~~Describe \code{expoints} here~~ }
  \item{xlab}{ ~~Describe \code{xlab} here~~ }
  \item{ylab}{ ~~Describe \code{ylab} here~~ }
  \item{lwd}{ ~~Describe \code{lwd} here~~ }
  \item{lty}{ ~~Describe \code{lty} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
%\references{  }
\author{ Robert Ferstl, Josef Hayden }
%\note{ }
\seealso{}
%\examples{}

\keyword{  }
