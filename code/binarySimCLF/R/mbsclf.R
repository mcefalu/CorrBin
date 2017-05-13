`mbsclf` <-
function(m, u, B, seed)
{
  if (!missing(seed))
  {
    set.seed(seed);
  }
  n <- nrow(B);
  y <- matrix(rep(-1, m*n), m);
  for (i in 1:m)
  {
    stage.i <- mbsclf1(u,B);
    if(!stage.i$succeed) {
        return(stage.i)
    }
    y[i,] <- t( stage.i$y ) ;
  }
  return( list( succeed=TRUE, y=y ) );
}

