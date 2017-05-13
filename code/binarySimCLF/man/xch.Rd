\name{xch}
\alias{xch}
\title{ Exchangeable Correlation Matrix }
\description{
    Constructs an \eqn{n} by \eqn{n} exchangeable correlation matrix.
}
\usage{
xch(n, rho)
}
\arguments{
  \item{n}{ Cluster size. }
  \item{rho}{ Exchangeable correlation parameter. }
}
\value{
    Returns an exchageable correlation matrix of dimension \code{n} by \code{n} with parameter \code{rho}.
}
\seealso{ See Also \code{\link{ar1}}, \code{\link{ma1}} }
\examples{
    #Example
    xch(5,0.17)
}
