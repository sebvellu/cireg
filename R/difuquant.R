difuquant <- function(
	prob, dpow = NULL, xnum = 0, step = 2, lgth = 1000, simu = 10000
) {
	dvls <- polydet(lgth, dpow)
	rslt <- matrix(NA_real_, simu, 2)
	#
	for (indx in 1:simu) {
		yvls <- cumsum(rnorm(lgth))
		if (xnum > 0) {
			xvls <- apply(matrix(rnorm(lgth * xnum), lgth), 2, cumsum)
			yvls <- qr.resid(qr(cbind(dvls, xvls)), yvls)
			dvls <- NULL
		} else {
			xvls <- NULL
		}
		objt <- difustat(yvls, dvls, xvls, step)
		rslt[indx, ] <- c(objt$cffs, objt$tsta)
	}
    #
	return(apply(rslt, 2, get_quantile, prob))
}
