`getBnds` <-
function(i, u, b)
{
  y <- (b > 0);
  nuMax <- u[i] + t( y - u[ 1:(i-1) ] ) %*% b ;
  y <- !y ;
  nuMin <- u[i] + t( y - u[1:(i-1)] ) %*% b ;
  return( list(nuMin = nuMin, nuMax=nuMax) );
}

