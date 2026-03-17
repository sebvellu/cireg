shinimcmprquant <- function(prob, zpow, zadd, lgth = 1000, simu = 10000) {
	if (is.matrix(zpow)) {
		xnum <- ncol(zpow) - 1
	} else {
		if (length(zpow) == 1) {
			xnum <- attr(zpow, "xnum")
		} else {
			xnum <- length(zpow) - 1
		}
	}
	rslt <- numeric(simu)
    for (indx in 1:simu) {
        yvls <- rnorm(lgth)
        xvls <- apply(matrix(rnorm(lgth * xnum), lgth), 2, cumsum)
        zvls <- getzvls(zpow, xvls)
		#
		s1zx <- cbind(apply(zvls, 2, cumsum), zadd, xvls)
		s1rs <- qr.resid(qr(s1zx), cumsum(yvls))
		rsds <- c(s1rs[1], diff(s1rs))
		#
        rslt[indx] <- kpssstatint(rsds, 1)
    }
	return(get_quantile(rslt, prob))
}
