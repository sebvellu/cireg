cmprreg <- function(zpow, xvls) {
	znum <- nrow(zpow)
	lgth <- nrow(xvls)
	#if (is.null(xvls)) {
	#	xnum <- 0
	#} else {
	#	xnum <- ncol(xvls)
	#}
	xnum <- ncol(xvls)
	base <- cbind(1:lgth, xvls)
	rslt <- matrix(0, lgth, znum)
	for (iidx in 1:znum) {
		temp  <- 1
		for (jidx in 1:(xnum + 1)) {
			temp <- temp * base[, jidx]^zpow[iidx, jidx]
		}
		rslt[, iidx] <- temp
	}
	return(rslt)
}
