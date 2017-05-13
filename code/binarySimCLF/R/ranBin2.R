`ranBin2` <-
function(nRep,u,psi,seed)
{
  if (!missing(seed))
  {
    set.seed(seed);
  }
  u12 <- solve2(u[1], u[2], psi);
  y <- matrix( rep(-1,2*nRep),nrow=nRep );
  y[,1] <- ifelse( runif(nRep)<=u[1], 1, 0 );
  y[,2] <- y[,1]*(runif(nRep) <= u12/u[1]) + (1-y[,1])*(runif(nRep) <= (u[2]-u12)/(1-u[1]));
  return(y)
}

