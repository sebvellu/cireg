ctimcmprquant <- function(prob, ynum, zpow, zadd, lgth = 1000, simu = 10000) { #zadd == NULL or 1 
	if (is.matrix(zpow)) {
		xnum <- ncol(zpow) - 1
	} else {
		if (length(zpow) == 1) {
			xnum <- attr(zpow, "xnum")
		} else {
			xnum <- length(zpow) - 1
		}
	}
	rslt <- matrix(0, simu, 3)
	#
    for (indx in 1:simu) {
        yvls <- matrix(stats::rnorm(lgth * ynum), lgth)
        xvls <- apply(matrix(stats::rnorm(lgth * xnum), lgth), 2, cumsum)
        zvls <- getzvls(zpow, xvls)
		#
		s1zx <- cbind(apply(zvls, 2, cumsum), zadd, xvls)
		s1yv <- apply(yvls, 2, cumsum)
		s1rs <- qr.resid(qr(s1zx), s1yv)
		#
		#vala <- crossprod(s1rs)
		#valb <- crossprod(apply(s1rs, 2, cumsum))
		#rslt[indx, 1] <- sum(diag(vala))/lgth^2
		#rslt[indx, 2] <- sum(diag(valb))/lgth^4
		#rslt[indx, 3] <- sum(diag(safesolve(vala, valb)))/lgth^2
		rslt[indx, ] <- ctstatint(s1rs, NULL)
	}
	#
	return(apply(rslt, 2, get_quantile, prob))
}