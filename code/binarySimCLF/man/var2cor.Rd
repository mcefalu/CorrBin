\name{var2cor}
\alias{var2cor}
\title{ Variance To Correlation }
\description{
    Converts a covariance matrix into a correlation matrix.
}
\usage{
var2cor(v)
}
\arguments{
  \item{v}{ A covariance matrix. }
}
\value{
  Returns a correlation matrix.
}

\seealso{ \code{\link{cor2var}} }
\examples{
    # Example
    v1 = c(25, 21, 17); v1
    v = toeplitz(v1); v
    var2cor(v)
}
