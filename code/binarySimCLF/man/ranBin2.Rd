\name{ranBin2}
\alias{ranBin2}
\title{ Bivariate Binary Data }
\description{
    Generates correlated binary vectors each of length 2.  The routine is based
    on the pairwise odds ratio.
}
\usage{
ranBin2(nRep, u, psi, seed)
}
\arguments{
  \item{nRep}{ Number of replications. One can think of \code{nRep} as the number
        of clusters with mean \code{u} and pairwise odds ratio \code{psi}. }
  \item{u}{ Mean vector of length 2 }
  \item{psi}{ Odds ratio. }
  \item{seed}{ Sets the seed. }
}
\value{
    Returns a response matrix of dimension \code{nRep} by 2.  Each row represents
    an independent cluster.
}
\seealso{ \code{\link{mbsclf1}}, \code{\link{mbsclf}}, \code{\link{ranXch}} }
\examples{
    # Examples
    u = c(0.25, 0.25);
    psi = 1.75;

    # Begin: checking whether the odds pairwise odds ratio is compatible with the
    # mean vector u.
    v = psi2var(psi,u);  v
    cr = var2cor(v);  cr
    cp = chkBinC(cr,u); cp$compat
    # End: checking compatibility.

    y = ranBin2(1000,u,psi)
    r = cor(y);r
    ybar = apply(y,2,mean);ybar
    psihat = cor2psi(r,ybar)$psi;psihat
}
