\name{var2psi}
\alias{var2psi}
\title{ From Covariance Matrix To Pairwise Odds Ratio }
\description{
    Computes the pairwise odds ratio vector, denoted by \eqn{\psi}, based on the
    covariance matrix \eqn{V} of a binary random vector with mean \eqn{\mu}.
}
\usage{
var2psi(v, mu)
}
\arguments{
  \item{v}{ Covariance matrix. }
  \item{mu}{ Mean vector. }
}
\value{
  Returns a list with two components:
  \item{psi}{The pairwise odds ratio vector if the routine executes without error.}
  \item{error}{A logical; if \code{TRUE} then an error has occurred. This usually means that
    division by zero has occured.}
}
\seealso{ \code{\link{cor2psi}} }
