\name{cor2psi}
\alias{cor2psi}
\title{ From Correlation To Pairwise Odds-Raio }
\description{
This function converts a correlation matrix associated with a mean
vector into a vector of pairwise odds ratio.  This is done in the following
order.  Let \eqn{\mu = (\mu_1, \mu_2, ..., \mu_n)}.  Then the odds ratio vector is of length
\eqn{n} choose 2 with elements \eqn{\psi_{12}, \psi_{13}, ..., \psi_{1n}, \psi_{23}, ..., \psi_{(n-1)n} }.
}
\usage{
cor2psi(r, mu)
}
\arguments{
  \item{r}{ Correlation matrix. }
  \item{mu}{ Mean vector. }
}
\value{
  Returns a list with two components:
  \item{psi}{The pairwise odds ratio vector if the routine executes without error.}
  \item{error}{A logical; if \code{TRUE} then an error has occurred. This usually means that
    division by zero has occured.}
}
\seealso{ \code{\link{var2psi}} }
\examples{
    # Example
    mu = c(0.2, 0.3, 0.25);
    R = xch(3,0.30);
    cor2psi(R,mu)
}