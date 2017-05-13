`ar1` <-
function(n, rho)
{
  if (n <= 0)
  {
    stop('n must be at least 1');
  }
  if (n <= 2)
  {
    return(xch(n,rho));
  }
  else
  {
    line1 <- c(1,rho^(1:(n-1)));
    return(toeplitz(line1));
  }
}

