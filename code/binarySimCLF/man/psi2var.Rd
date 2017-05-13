\name{psi2var}
\alias{psi2var}
\title{ From Odds Ratio To Covariance  }
\description{
    Takes a vector of pairwise odds ratio (in lexicographic order) and gives
    an equivalent covariance matrix.
}
\usage{
psi2var(psi, mu)
}
\arguments{
  \item{psi}{ Vector of pairwise odds ratio. }
  \item{mu}{ Mean vector. }
}
\details{
    If we have a cluster of size 4, then the \code{psi} vector must have elements of
    the form \eqn{[\psi_{12}, \psi_{13}, \psi_{14}, \psi_{23}, .., \psi_{34}]^T}.
}
\value{
    Returns an equivalent covariance matrix.
}
\examples{
    # Example
    r = toeplitz(c(1, 0.45, 0.37));
    mu = c(0.25, 0.32, 0.32); mu

    # Checks to see whether the specified mean and correlation matrix are compatible.
    # The correlation parameters are restricted by functions of the pairwise odds
    # ratio.
    cp = chkBinC(r,mu);cp

    # If r and mu are compatible, compute the covariance matrix.
    if(cp$compat) { v = cor2var(r,mu);v }

    # Computes the pairwise odds ratio.
    or <- var2psi(v,mu); or

    # Converts from the pairwise odds ratio to the covariance matrix.
    psi2var(or$psi, mu)
}
