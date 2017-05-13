\name{rhoRange}
\alias{rhoRange}
\title{ Range Of Correlation }
\description{
    Given a mean vector \eqn{\mu} of length \eqn{n}, the function computes the lower and upper
    bound for an exchangeable correlation parameter imposed by the mean and the cluster
    size.
}
\usage{
rhoRange(u)
}
\arguments{
  \item{u}{ Mean vector. }
}

\value{
  Returns a list with two components.
  \item{rhomin }{The smallest value that the correlation parameter is permitted to take.}
  \item{rhomax }{The largest value that the correlation parameter is permitted to take.}
}
\seealso{\code{\link{ranXch}} }
\examples{
    # Examples
    u=c(0.25, 0.12, 0.27, 0.76, 0.95);u
    u1 = c(0.35, 0.42);u1
    rhoRange(u)
    rhoRange(u1)
    rhoRange(0.25)
}

