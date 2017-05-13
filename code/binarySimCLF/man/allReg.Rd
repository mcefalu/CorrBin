\name{allReg}
\alias{allReg}
\title{ An Intermediate Matrix }
\description{
	This function computes an INTERMEDIATE matrix that is to be used by other routines.
	Details of this intermediate matrix is found in Qaqish (2003).  This intermediate
	is needed for the actual simulation.
}
\usage{
allReg(V)
}

\arguments{
  \item{V}{ The covariance matrix cov(Y). }
}

\value{
  Returns a matrix that is some function of V.
}

\examples{
    #Examples

    # correlation matrix
    r = cbind( c(1, 0.1, 0.2), c(0.1, 1, 0.3), c(0.2, 0.3, 1) );

    # mean vector
    mu = c(0.4, 0.5, 0.6);

    # covariance matrix
    v = cor2var(r, mu) ;

    # Intermediate matrix
    B = allReg(v);
}
