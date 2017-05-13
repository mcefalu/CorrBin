\name{getBnds}
\alias{getBnds}
\title{ Bounds Of Conditional Mean }
\description{
    Let \eqn{\lambda_i = E(Y_i|Y_1, ..., Y_{i-1})}.  For a given \eqn{i}, the function
    computes the bounds for \eqn{\lambda_i}.  There is no reason that the user should
    use this function directly.  It is called by the function \code{blrchk1}
    which is used to check for compatibility with the CLF.
}
\usage{
getBnds(i, u, b)
}
\arguments{
  \item{i}{ Indexes the \eqn{i}-th response. }
  \item{u}{ Mean vector. }
  \item{b}{ Vector of regression coefficients of length \eqn{i-1}. }
}
\details{
    The vector \code{b} is obtained from the columns of B = \code{allReg(V)}.  See expressions
    (3) and Lemma 1 in Qaqish(2003).
}
\value{
    Returns a list with two components:
  \item{nuMin}{Lower bound for \eqn{\lambda_i}.}
  \item{nuMax}{Upper bound for \eqn{\lambda_i}.}
}
\seealso{\code{\link{blrchk1}} }

