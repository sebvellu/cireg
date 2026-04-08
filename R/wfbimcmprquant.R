wfbimcmprquant <- function(
	prob, ynum, zpow, hypo, krnl, bfrc, zadd, lgth = 1000, simu = 10000
) { #zadd == NULL or 1 
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
	band <- ceiling(bfrc * lgth)
	#
    for (indx in 1:simu) {
        yvls <- matrix(stats::rnorm(lgth * ynum), lgth)
        xvls <- apply(matrix(stats::rnorm(lgth * xnum), lgth), 2, cumsum)
        zvls <- getzvls(zpow, xvls)
		#
		# s1zx <- cbind(apply(zvls, 2, cumsum), zadd, xvls)
		# s1yv <- apply(yvls, 2, cumsum)
		# #
		# orth <- apply(s1zx[lgth:1, , drop = FALSE], 2, cumsum)
		# orth <- apply(orth[lgth:1, , drop = FALSE], 2, cumsum)
		# # equivalently #
		# #s2zx <- apply(s1zx, 2, cumsum)
		# #cvls <- sweep(-s2zx[1:(lgth - 1), , drop = FALSE], 2, colSums(s1zx), "+")
		# #cvls <- rbind(colSums(s1zx), cvls)
		# #orth <- apply(cvls, 2, cumsum)
		# s1zo <- cbind(s1zx, orth)
		# rsdm <- qr.resid(qr(s1zo), s1yv)
		rsdm <- imlsfit(yvls, zvls, xvls, zadd, NULL, NULL, NULL)$rsdm
		lrvr <- lrvar::lrvar(rbind(rsdm[1, ], diff(rsdm)), krnl, band)$longvar
		nvls <- stats::rnorm(hypo * ynum)
		lrvi <- helperkit::safesolve(lrvr)
		stat <- sum(nvls * ((lrvi %x% diag(1, hypo)) %*% nvls))
		#
		rslt[indx] <- stat
	}
	#
	return(helperkit::get_quantile(rslt, prob))
}