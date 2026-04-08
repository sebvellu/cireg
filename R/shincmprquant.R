shincmprquant <- function(prob, zpow, lgth = 1000, simu = 10000) {
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
        yvls <- stats::rnorm(lgth)
        xvls <- apply(matrix(stats::rnorm(lgth * xnum), lgth, xnum), 2, cumsum)
        rsds <- qr.resid(qr(getzvls(zpow, xvls)), yvls)
        rslt[indx] <- kpssstatint(rsds, 1)
    }
	return(helperkit::get_quantile(rslt, prob))
}
