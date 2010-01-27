\name{grad_objfct_sv}
\Rdversion{1.1}
\alias{grad_objfct_sv}

\title{Gradient of the Svensson Loss Function for Yields
}
\description{
Calculates the gradient of the Svensson Loss Function for Yields
}
\usage{
grad_objfct_sv(beta, m, y)
}

\arguments{
  \item{beta}{Parameter of the Svensson spot rate function (for
  details see \code{\link{spr_sv}}). 
}
  \item{m}{maturity vector
}
  \item{y}{yield vector
}
}


\seealso{ \code{\link{objfct_sv}},\code{\link{spr_sv}}
}