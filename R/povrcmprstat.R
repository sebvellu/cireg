povrcmprstat <- function(
    yvls, zpow, xvls, krnl = "ba", band = "and", long = NULL
) {
	if (is.matrix(zpow)) {
		didx <- apply(zpow[, -1, drop = FALSE], 1, function(x) {
            return(all(x == 0))
        })
	} else {
		if (length(zpow) == 1) {
            didx <- c(rep(TRUE, zpow[1] + 1), rep(FALSE, ncol(xvls)))
		} else {
            didx <- c(rep(TRUE, zpow[1] + 1), rep(FALSE, sum(zpow[-1])))
		}
	}
	zvls <- getzvls(zpow, xvls)
	dvls <- zvls[, didx, drop = FALSE]
	xaug <- zvls[, !didx, drop = FALSE]
    #
    return(drop(povrstat(yvls, xvls, dvls, xaug, krnl, band, long)))
}
