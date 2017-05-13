`var2cor` <-
function(v)
{
  p = dim(v)[1];
  d = 1/(sqrt(diag(v)));
  return( d * v * rep(d,each=p));
}

