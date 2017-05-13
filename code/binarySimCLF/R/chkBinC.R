`chkBinC` <-
function(r, mu)
{
  n <- length(mu);  # use LENGTH() for vectors instead of NROW()
  v <- cor2var(r, mu);
  for (i in 1:(n-1) )
  {
    mui <- mu[i];
    for (j in (i+1):n)
    {
      muj <- mu[j];
      muij <- mui*muj + v[i,j];
      ok <- ( muij <= min(mui, muj) ) & ( muij >= max(0, mui + muj - 1) );
      if (!ok)
      {
        return(list(compat=ok, locOfFail = paste("i = ",i,", j = ", j)));
      }
    }
  }
  return( list(compat = ok, locFail = paste("Compatible:  no failure") ) );
}

