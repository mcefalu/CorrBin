`rhoRange` <-
function(u)
{
  n <- length(u);
  if (n == 1)
  {
    return(list(rhomin=-1, rhomax=1)) ;
  }
  u1 <- min(u);
  un <- max(u);
  if (n==2)
  {
    u2 <- un;
    un1 <- u1;
  }
  else
  {
    # Finds the second smallest (stored in u2)
    i <- (u == u1);
    u[i] <- un;
    u2 <- min(u);
    u[i] <- u1;
    # Finds the second largest (stored in un1)
    i <- (u == un);
    u[i] <- u1 ;
    un1 <- max(u) ;
    u[i] <- un ;
  }
  psi1 <- sqrt(u1/(1-u1));
  psi2 <- sqrt(u2/(1-u2));
  psin <- sqrt(un/(1-un));
  psin1 <- sqrt(un1/(1-un1));

  rhomin = -min(psi1*psi2, 1/(psin1*psin), 1/(n-1));
  rhomax = psi1/psin;
  return(list(rhomin=rhomin, rhomax=rhomax))
}

