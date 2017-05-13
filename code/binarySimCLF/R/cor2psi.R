`cor2psi` <-
function(r,mu)
{
  n <- nrow(r)
  psi <- rep( -1, choose(n,2) )
  v <- cor2var(r,mu)
  k <- 1
  for (i in 1:(n-1))
  {
    mui <- mu[i]
    for (j in (i+1):n)
    {
      muj <- mu[j]
      muij <- mui*muj + v[i,j]
      psi[k] <- muij * (1 - mui - muj + muij)/( (mui-muij)*(muj - muij) )
      k <- k+1
    }
  }
  return(list(psi = psi, error = any(!is.finite(psi))))
}

