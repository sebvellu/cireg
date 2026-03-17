cmprdzdx <- function(xvls, zpow) { #parm, xvls, zvls, xexp, zexp) { #xvls zpow
	lgth <- nrow(xvls)
	base <- cbind(1:lgth, xvls)
	xfct <- rep(1, ncol(xvls))
	rslt <- array(dim = c(nrow(zpow), ncol(xvls), lgth))
	for (i in 1:nrow(zpow)) {
		xpow <- zpow[i, -1]
		tpow <- zpow[i, 1]
		for (j in 1:ncol(xvls)) {
			fctx <- replace(xfct, j, xpow[j])
			powx <- replace(xpow, j, ifelse(xpow[j] < 1, 0, xpow[j] - 1))
			fctz <- rep(c(1, fctx), each = lgth)
			powz <- rep(c(tpow, powx), each = lgth)
			rslt[i, j, ] <- apply(fctz * base^powz, 1, prod)
		}
	}
	return(rslt)
}