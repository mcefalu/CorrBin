\name{blrchk}
\alias{blrchk}
\title{ CLF Compatibility 2 }
\description{
    Checks for CLF compatibility.  Unlike \code{blrchk1()}, this function
    works on the covariance matrix \eqn{V}.  It calls \code{blrchk1()}.  It is
    RECOMMENDED that the user use this function or \code{blrchk1()} first before
    engaging in the actual simulation.
}
\usage{
blrchk(u, V)
}

\arguments{
  \item{u}{ Mean vector }
  \item{V}{ Covariance matrix }
}

\value{
    Returns a logical indicating CLF compatibility; \code{TRUE} denoting
    CLF compatibility and \code{FALSE} otherwise.
}
\seealso{ \code{\link{blrchk1}}, \code{\link{binarySimCLF}}. }
\examples{
    # Example
    # Checks CLF compatibility.
    R = xch(3,0.9);
    mu = c(0.4,0.5,0.6);

    V = cor2var(R,mu); V
    clf.compat = blrchk(mu,V); clf.compat
}
