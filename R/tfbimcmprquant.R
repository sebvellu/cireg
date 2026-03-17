tfbimcmprquant <- function(
	prob, zpow, krnl, bfrc, zadd, lgth = 1000, simu = 10000
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
        yvls <- rnorm(lgth)
        xvls <- apply(matrix(rnorm(lgth * xnum), lgth), 2, cumsum)
        zvls <- getzvls(zpow, xvls)
		#
		# s1zx <- cbind(apply(zvls, 2, cumsum), zadd, xvls)
		# s1yv <- cumsum(yvls)
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
		lrvr <- lrvar(c(rsdm[1], diff(rsdm)), krnl, band)$longvar
		#
		rslt[indx] <-  rnorm(1)/sqrt(lrvr)
	}
	#
	return(get_quantile(rslt, prob))
}