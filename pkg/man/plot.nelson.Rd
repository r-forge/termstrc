\name{plot.nelson}
\alias{plot.nelson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.nelson(x, matrange = c(min(mapply(function(i) min(x$y[[i]][, 1]), seq(x$n_group))), max(mapply(function(i) max(x$y[[i]][, 1]), seq(x$n_group)))), multiple = FALSE, expoints = unlist(x$expoints), ctype = "spot", errors = "price", lwd = 2, lty = 1, type = "l", inset = c(0.8, 0.1), ask = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{matrange}{ ~~Describe \code{matrange} here~~ }
  \item{multiple}{ ~~Describe \code{multiple} here~~ }
  \item{expoints}{ ~~Describe \code{expoints} here~~ }
  \item{ctype}{ ~~Describe \code{ctype} here~~ }
  \item{errors}{ ~~Describe \code{errors} here~~ }
  \item{lwd}{ ~~Describe \code{lwd} here~~ }
  \item{lty}{ ~~Describe \code{lty} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{inset}{ ~~Describe \code{inset} here~~ }
  \item{ask}{ ~~Describe \code{ask} here~~ }
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

