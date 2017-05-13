`allReg` <-
function(V)
{
  n <- nrow(V);
  B <- V;
  for (i in 2:n)
  {
    i1 <- i-1;
    gi <- V[1:i1, 1:i1];
    si <- V[1:i1, i];
    Bi <- solve(gi, si);
    B[1:i1, i] <- Bi ;
  }
  return(B);
}

