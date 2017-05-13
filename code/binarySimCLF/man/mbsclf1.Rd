\name{mbsclf1}
\alias{mbsclf1}
\title{ One Correlated Binary Vector }
\description{
    This function generates a single correlated binary vector.  It is called by
    \code{mbsclf()} to generate multiple vectors of correlated binary data.
}
\usage{
mbsclf1(u, B)
}
\arguments{
  \item{u}{ Mean vector. }
  \item{B}{ Intermediate vector obtained by \code{B = allReg(V)} .}
}
\value{
    Returns a list with two components.
    \item{succeed }{Indicates whether the routine succeeds or fails.}
    \item{y }{The response.  If \code{succeed} is \code{FALSE}, then \code{y}
              is \code{NULL}.  If \code{succeed} is \code{TRUE}, then \code{y}
              is a binary outcome vector.}
}
\details{
Failure means that the condtional mean is out of range.  Using the
notation of Qaqish (2003), it means that \eqn{\lambda_i} is not in the
interval \eqn{[0,1]}.  This is one of the criteria for CLF compatibility.
Note also that even if the code succeeds, this does not mean that
we have CLF compatibility.  To ensure CLF compatibility, we have
to make sure that the largest possible value of \eqn{\lambda_i} is less than or
equal to 1 and the smallest possible value of \eqn{\lambda_i} is greater than or
equal to 0.  These concepts are explained in Qaqish (2003).
}
\references{Qaqish, B. F., A family of multivariate binary distributions for
simulating correlated binary variables with specified marginal means and
correlations. \emph{Biometrika} \bold{92}:455-463, 2003.
}
\seealso{ \code{\link{mbsclf}}, \code{\link{allReg}}, \code{\link{binarySimCLF}}. }
\examples{
    #Examples
    R = cbind( c(1, 0.1, 0.2), c(0.1, 1, 0.3), c(0.2, 0.3, 1) );
    mu = c(0.4, 0.5, 0.6);
    V = cor2var(R, mu) ;
    B = allReg(V);
    clf.compat = blrchk(mu,B);
    if (clf.compat)
    {
        y = mbsclf1(mu, B);
        print("CLF compatible");
    } else { print("Not CLF compatible"); }
}
