\name{mardia}
\alias{mardia}
\title{ Quadratic Equation Solver }
\description{
    Consider the quadratic equation: \eqn{a x^2 + bx + c = 0}.  This function
    returns the appropriate zero of this equation based on Mardia's formula.
    The user should not use this function directly.  It is called by \code{solve2()}.
}
\usage{
mardia(a, b, c)
}
\arguments{
  \item{a}{ Coefficient of \eqn{x^2}.  Should be nonzero. }
  \item{b}{ Coefficient of \eqn{x}. }
  \item{c}{ Constant. }
}
\value{
    Returns the appropriate root of the quadratic equation.  Note that only one
    root is returned.
}
\examples{
    # Example
    mardia(1,4,3)
    mardia(1,-4,3)
}
