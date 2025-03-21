inla.ar.pacf2phi<-function (pac) 
{
  p <- length(pac)
  stopifnot(p > 0)
  phi <- work <- pac
  if (p > 1L) {
    for (j in 1L:(p - 1L)) {
      a <- phi[j + 1L]
      phi[1L:j] <- work[1L:j] <- work[1L:j] - a * phi[j:1L]
    }
  }
  return(phi)
}