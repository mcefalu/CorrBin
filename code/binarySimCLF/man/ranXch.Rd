\name{ranXch}
\alias{ranXch}
\title{ Correlated Binary Data With Exchangeable Correlation }
\description{
  Generates correlated binary data with an exchangeable correlation.
}
\usage{
ranXch(u, alpha, seed)
}
\arguments{
  \item{u}{ The mean vector. }
  \item{alpha}{ The exchangeable correlation parameter. }
  \item{seed}{ Sets the seed. }
}
\value{
    The function returns a list with two components.
    \item{succeed}{A logical indicating whether the routine succeeded.  \code{TRUE} if
        it succeeded.  \code{FALSE} if it did not.  \code{FALSE} indicates that it is NOT CLF
        compatible.}
    \item{y}{Outcome: If \code{succeed} is \code{FALSE} then \code{y} is \code{NULL}.
        If \code{succeed} is \code{TRUE} then \code{y} is the outcome vector.}
}
\seealso{ \code{\link{binarySimCLF}}. }
\examples{
    #Example
    u = c(0.25, 0.24, 0.26, 0.35); u
    ranXch(u,0.25)
}
