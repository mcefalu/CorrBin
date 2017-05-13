\name{chkBinC}
\alias{chkBinC}
\title{ Compatibility of the Mean and Correlation }
\description{
    Checks whether the pairwise correlation is compatible with the mean vector.
    The correlation parameter has to satisfy some natural restrictions and this
    function checks whether this is the case.  Note that this is not the same
    thing as CLF compatibility.
}
\usage{
chkBinC(r, mu)
}
\arguments{
  \item{r}{ Correlation matrix. }
  \item{mu}{ Mean vector. }
}

\value{
    Returns a list with two components.
  \item{compat}{A logical such that if \code{TRUE}, then there is compatibility.}
  \item{locFail}{If not compatible, this is a string giving the location of the failure.
    If compatible, this is a string telling the user that the mean and the
    correlation is compatible. }
}

\examples{
    # Examples

    R = cbind( c(1, 0.1, 0.2), c(0.1, 1, 0.3), c(0.2, 0.3, 1) )
    mu = c(0.4, 0.5, 0.6)
    v = cor2var(R,mu)
    natural.compat = chkBinC(R,mu); natural.compat
}
