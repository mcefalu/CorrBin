`ranXch` <-
function(u, alpha, seed)
{
  n <- length(u)
  std <- sqrt(u*(1-u))
  a <- 1/std
  d <- 1:n - 2
  d <- alpha / (1 + alpha*d)
  dstd <- d * std
  y <- rep(-1, n)
  if (!missing(seed))
  {
    set.seed(seed);
  }
  y[1] <- ifelse( runif(1) <= u[1],1,0 );
  r <- rep(0,n);
  r[1] <- y[1] - u[1];
  lam <- u;
  for (i in 2:n)
  {
    j <- 1:(i-1);  # Note:  j is a vector used to subset other vectors.
    bi <- dstd[i] * a[j];
    lam[i] <- u[i] + t(r[j])%*%bi  ;
    y[i] <- ifelse( runif(1) <= lam[i], 1, 0 );
    r[i] <- y[i] - u[i] ;
  }

  if ((min(lam) < 0) | (max(lam) > 1))
  {
    # print(paste("Conditional means outside [0,1]) ",lam));
    return( list(succeed=FALSE, y=NULL) );
  }
  return( list(succeed=TRUE, y=y) )
}

