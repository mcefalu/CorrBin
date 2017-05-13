\name{ma1}
\alias{ma1}
\title{ Moving Average Of Order 1 }
\description{
    Constructs an \eqn{n} by \eqn{n} MA(1) correlation matrix.
}
\usage{
ma1(n, rho)
}
\arguments{
  \item{n}{ Dimension of correlation matrix. }
  \item{rho}{ Correlation parameter. }
}
\value{
  Returns an \code{n} by \code{n} MA(1) correlation matrix.
}
\seealso{ \code{\link{ar1}}, \code{\link{xch}}}
\examples{
    #Example
    ma1(5,0.56)
}
