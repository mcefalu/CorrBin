\name{solve2}
\alias{solve2}
\title{ Pairwise Moments }
\description{
    Given two marginal means mui, muj and a pairwise odds ratio psi, the function
    computes the moment \eqn{E[Y_i Y_j] := \mu_{ij}}.  This function is called by
    \code{ranBin2} and \code{psi2var}.  There is no reason for the user to use this
    function in any direct way.
}
\usage{
solve2(mui, muj, psi)
}
\arguments{
  \item{mui}{ Mean of Bernoulli \eqn{Y_i} }
  \item{muj}{ Mean of Bernoulli \eqn{Y_j} }
  \item{psi}{ Pairwise odds ratio of \eqn{Y_i} and \eqn{Y_i} }
}
\value{
    Returns \eqn{E[Y_i Y_j] := \mu_{ij}}.
}
\seealso{ \code{\link{mardia}}, \code{\link{ranBin2}}, \code{\link{psi2var}} }
\examples{
    #Example
    mu1=0.35;
    mu2 = 0.35;
    psi = 2.5 ;
    solve2(mui=mu1, muj=mu2, psi=2.5)
}
