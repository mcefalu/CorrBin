\name{blrchk1}
\alias{blrchk1}
\title{ CLF Compatibility 1 }
\description{
    Checks for CLF compatibility given an intermediate matrix \eqn{B} which is obtained
    by \code{B = allReg(V)} where \eqn{V} is the covariance matrix.  If
    the intermediate matrix \eqn{B} is already computed, then it is best
    to use this routine.  If \eqn{B} is not precomputed then use \code{blrchk}
    which requires you to give it the covariance matrix \eqn{V}.
}
\usage{
    blrchk1(u, B)
}

\arguments{
  \item{u}{ Mean vector. }
  \item{B}{ An intermediate matrix. }
}

\value{
    Returns a logical indicating CLF compatibility: \code{TRUE} for compatible and
    \code{FALSE} for not compatible.
}

\seealso{ \code{\link{blrchk}}, \code{\link{binarySimCLF}}. }
\examples{
# Checks CLF compatibility.
    R = xch(3,0.9);
    mu = c(0.4,0.5,0.6);

    V = cor2var(R,mu); V
    B=allReg(V);
    clf.compat = blrchk1(mu,B); clf.compat
}
