`blrchk1` <-
function(u, B)
{
  n <- nrow(B) ;
  rc <- TRUE ;
  i <- 2;
  
  # do-while loop:
  while (rc  & ( i < (n+1) ) )
  {
    bounds <- getBnds( i, u, B[1:(i-1), i] ) ;
    rc <- ( (bounds$nuMax <= 1.0) & (bounds$nuMin >= 0.0) );
    i <- i+1;
  }

  return(as.vector(rc));
}

