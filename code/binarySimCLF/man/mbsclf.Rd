\name{mbsclf}
\alias{mbsclf}
\title{ Multiple Correlated Binary Vectors }
\description{
    Unlike \code{mbsclf1}, this routine generates multiple vectors where each
    row of the outcome matrix corresponds to a cluster.

    Please be advised that this function DOES NOT check for CLF compatibility.
    The onus is on the user to do so.  His or her program MUST check for CLF
    compatibility, otherwise the returned value is nonsense.  The routine for
    checking CLF compatibility is \code{blrchk}.
}
\usage{
mbsclf(m, u, B, seed)
}

\arguments{
  \item{m}{ The number of vectors of correlated binary outcome.  All have the
            same mean and correlation.  Thus, the generated vectors of responses are
            iid.  From a clustered data perspective, think of \code{m} as the number
            of clusters.}
  \item{u}{ Mean vector. }
  \item{B}{ Intermediate matrix obtained by \code{B = allReg(V)}. }
  \item{seed}{ Sets the seed. }
}
\value{
    Returns a list with two components.
    \item{succeed }{Indicates whether the routine succeeds or fails.}
    \item{y }{The response.  If \code{succeed} is \code{FALSE}, then \code{y}
              is \code{NA}.  If \code{succeed} is \code{TRUE}, then \code{y}
              is a binary outcome matrix.  The number of rows of this matrix
              is \code{m}.}
}
\details{
As we said in \code{mbsclf1}, CLF compatibility requires that \eqn{\lambda_i} is in
\eqn{[0,1]} and also that the largest and smallest possible value of \eqn{\lambda_i} is in
range.  As we said in the description section, the latter compatibility condition
is not checked by this function.  It is checked by the function \code{blrchk}.
If \code{blrchk} returns a \code{FALSE}, then this routine should not be used.  If the
user does not check, then for some \eqn{i} out of the \code{m} vectors, \eqn{\lambda_i} may
be in range and for some other \eqn{i'} it may not.  This constitutes a failure of
this function.  The code is allowed to run until a failure is encountered at which
point it stops and returns a \code{NULL} value for the outcome matrix.  Note that if the
smallest and largest possible \eqn{\lambda_i} is out of range, it is possible for the routine
to succeed.  For large \code{m}, this is very unlikely.  Even if it did succeed, the
output is still nonsense.

The take home message is check for CLF compatibility using \code{blrchk} before
allowing either \code{mbsclf} or \code{mbsclf1} to run.  This ensures CLF
compatibility.
}
\references{Qaqish, B. F., A family of multivariate binary distributions for
simulating correlated binary variables with specified marginal means and
correlations. \emph{Biometrika} \bold{92}:455-463, 2003.
}
\seealso{ \code{\link{mbsclf1}}, \code{\link{binarySimCLF}}. }
\examples{
    # Examples
    R = cbind( c(1, 0.1, 0.2), c(0.1, 1, 0.3), c(0.2, 0.3, 1) );
    mu = c(0.4, 0.5, 0.6);
    V = cor2var(R, mu) ; V
    # Check CLF compatibility.
    clf.compat = blrchk(mu,V);
    
    if (clf.compat){
        print("CLF compatible");
        B = allReg(V);
        y = mbsclf(200,mu,B, seed=568381)
    }
    # sample mean
    colMeans(y$y);
    cor(y$y);
}
