inla.ar.phi2pacf<-function(phi){
  p <- length(phi)
  stopifnot(p > 0)
  work <- pac <- phi
  if (p > 1L) {
    for (j in (p - 1L):1L) {
      a <- pac[j + 1L]
      pac[1L:j] <- work[1L:j] <- (pac[1L:j] + a * pac[j:1L])/(1 - 
                                                                a^2)
    }
  }
  return(pac)
}