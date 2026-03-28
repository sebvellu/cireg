povrcmprquant <- function(prob, zpow, lgth = 1000, simu = 10000) {
	if (is.matrix(zpow)) {
		xnum <- ncol(zpow) - 1
		didx <- apply(zpow[, -1, drop = FALSE], 1, function(x) {
            return(all(x == 0))
        })
	} else {
		if (length(zpow) == 1) {
			xnum <- ncol(xvls)
            didx <- c(rep(TRUE, zpow[1] + 1), rep(FALSE, xnum))
		} else {
			xnum <- length(zpow) - 1
            didx <- c(rep(TRUE, zpow[1] + 1), rep(FALSE, sum(zpow[-1])))
		}
	}
	rslt <- numeric(simu)
	for (indx in 1:simu) {
		yvls <- cumsum(stats::rnorm(lgth))
		xvls <- apply(matrix(stats::rnorm(lgth * xnum), lgth), 2, cumsum)
		zvls <- getzvls(zpow, xvls)
		dvls <- zvls[, didx, drop = FALSE]
		xaug <- zvls[, !didx, drop = FALSE]
		rslt[indx] <- povrstat(yvls, xvls, dvls, xaug, NULL, NULL, 1)
	}
	return(get_quantile(rslt, prob))
}