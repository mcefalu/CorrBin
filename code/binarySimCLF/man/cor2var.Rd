\name{cor2var}
\alias{cor2var}
\title{ From Correlation To Covariance }
\description{
    Converts a correlation matrix to a covariance matrix.
}
\usage{
cor2var(r, mu)
}
\arguments{
  \item{r}{ Correlation matrix. }
  \item{mu}{Mean vector. }
}
\value{
    Returns a covariance matrix.
}
\examples{
    # Examples
    r = toeplitz(c(1, 0.45, 0.37)); r
    mu = c(0.25, 0.32, 0.79); mu
    cor2var(r,mu)
    var2cor(cor2var(r,mu))
}