cprtocmprpow <- function(dpow, xpow) {
  lgth <- length(xpow) + 1
  rslt <- matrix(0, nrow = sum(xpow) + dpow + 1, ncol = lgth)
  posi <- dpow + 2
  if (dpow > -1) {
    rslt[1:(dpow + 1), 1] <- 0:dpow
  } 
  for (iidx in seq_along(xpow)) {
    xidx <- xpow[iidx]
    rslt[posi:(posi + xidx - 1), iidx + 1] <- seq_len(xidx)
    posi <- posi + xidx
  }
  return(rslt)
}
