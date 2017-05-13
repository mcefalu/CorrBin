`psi2var` <-
function(psi, mu)
{
  n <- length(mu); # use LENGTH() for vectors instead of NROW()
  v <- diag(mu*(1-mu));
  k <- 1;
  for (i in 1:(n-1))
    {
      mui <- mu[i];
      for (j in (i+1):n)
      {
        v[i,j] <- solve2(mui, mu[j], psi[k]) - mui*mu[j];
        v[j,i] <- v[i,j];
        k <- k+1;
      }
    }
    return(v);
}

