`cor2var` <-
function(r, mu)
{
  p <- dim(r)[1];
  d <- sqrt(mu*(1-mu));
  V <- d * r * rep(d,each=p);
  rownames(V) <- NULL;
  return(V);
}

