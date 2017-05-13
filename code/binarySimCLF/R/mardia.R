`mardia` <-
function(a,b,c)
{
  if (a == 0 )
  {
    return(-c/b);

  }
  k <- ifelse( b > 0, 1, ifelse(b < 0, -1, 0) );
  p <- -0.5 * (b + k*sqrt(b^2-4*a*c));
  r1 <- p/a;
  r2 <- c/p;

  r <- ifelse( r2 > 0, r2, r1 );
  return(r);
}

