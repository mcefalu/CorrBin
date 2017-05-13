\name{ar1}
\alias{ar1}
\title{ Construct An AR(1) Correlation Matrix }
\description{
  This function constructs an AR(1) correlation matrix.
}
\usage{
ar1(n, rho)
}

\arguments{
  \item{n}{ The dimension of the correlation matrix. This usually corresponds
            to cluster size.}
  \item{rho}{ The correlation parameter. }
}

\value{
    Returns an AR(1) correlation matrix satisfying \code{n} and \code{rho}.
}

\seealso{ \code{\link{ma1}}, \code{\link{xch}}  }
\examples{
    # Example
    ar1(5,0.75)
    ar1(3,0.25)
}

