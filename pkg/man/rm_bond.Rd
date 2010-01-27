\name{rm_bond}
\Rdversion{1.1}
\alias{rm_bond}

\title{Remove Bonds from a Static Data set }
\description{
 Specified bonds and their associated data are removed from a bond data set.
}
\usage{
rm_bond(bonddata, ISIN, group)
}

\arguments{
  \item{bonddata}{bond data set. }
  \item{ISIN}{the ISIN numbers of the bonds to remove. }
  \item{group}{the group where the bonds to be removed belong to.}
}

\value{The function returns the new bond data set.}

\examples{
data(govbonds) 
newgovbonds <- rm_bond(govbonds,"DE0001135150","GERMANY")
}