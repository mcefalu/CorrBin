`solve2` <-
function(mui, muj, psi)
{
  if (psi == 1)
  {
    return(mui*muj);
  }
  else if (psi != 1)
  {
    a <- 1 - psi;
    b <- 1 - a*(mui + muj);
    c <- -psi*(mui * muj);
    muij <- mardia(a,b,c);
  }
  return(muij);
}

