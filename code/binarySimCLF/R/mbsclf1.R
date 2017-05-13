`mbsclf1` <-
function(u, B)
{
  n <- length(u);
  y <- rep(-1, n);
  y[1] <- ( runif(1) <= u[1] );
  r <- t(y - u);
  for (i in 2:n)
  {
    i1 <- i-1;
    ci <- u[i] + (r[1, 1:i1] %*% B[1:i1,i]);
    if ((ci < 0) || (ci > 1))
    {
      # print(paste("ERROR: Conditional mean outside of [0,1] in routine mbsclf1", ci ) );
      return( list(succeed=FALSE, y=NULL) );
      # Failure to succeed means that the conditional mean is outside of the
      # interval [0,1].
    }
    y[i] <- ( runif(1) <= ci );
    r[i] <- ( y[i] - u[i] ) ;
  }
  return( list(succeed=TRUE, y=y) );
}

