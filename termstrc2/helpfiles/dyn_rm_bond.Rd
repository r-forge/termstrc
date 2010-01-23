\name{dyn_rm_bond}
\Rdversion{1.1}
\alias{dyn_rm_bond}

\title{Remove Bonds from a Dynamic Data set}

\description{Specified bonds and their associated data are removed from a bond data set.}

\usage{
dyn_rm_bond(dynbonddata, group, ISIN)
}


\arguments{
  \item{dynbonddata}{dynamic bond data set.}
  \item{group}{the group where the bonds to be removed belong to.}
  \item{ISIN}{the ISIN numbers of the bonds to remove.}
}

\value{The function returns the new bond data set.}


\examples{
data(GermanBonds)
newdynbonddata <- dyn_rm_bond(datadyncouponbonds,"GERMANY","DE0001135150")
}

