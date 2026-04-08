phzfbquant <- function(
	prob, dpow = NULL, xnum = 0, step = 2, krnl = "ba", bfrc = 0.02, 
	lgth = 1000, simu = 10000
) {
	band <- ceiling(bfrc * lgth)
	#
	dvls <- helperkit::polydet(lgth, dpow)
	rslt <- matrix(NA_real_, simu, 2) #matrix(0, simu, 4)
	#
	for (indx in 1:simu) {
		yvls <- cumsum(stats::rnorm(lgth))
		if (xnum > 0) {
			xvls <- apply(matrix(stats::rnorm(lgth * xnum), lgth), 2, cumsum)
			yvls <- qr.resid(qr(cbind(dvls, xvls)), yvls)
			dvls <- NULL
		} else {
			xvls <- NULL
		}
		objt <- difustat(yvls, dvls, xvls, step)
		rslt[indx, ] <- phzfbstat(objt, krnl, band)[1:2]
	}
    #
	return(apply(rslt, 2, helperkit::get_quantile, prob))
}
