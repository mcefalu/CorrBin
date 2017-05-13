`xch` <-
function(n,rho)
{
  if (n <= 0)
  {
    stop('n must be at least 1');
  }
  else if (n==1)
  {
    r <- 1;
  }
  else if (n > 1)
  {
    c1 = rep(rho,n-1);
    r = toeplitz(c(1,c1)) ;
  }
  return(r);
}

