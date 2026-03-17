povrstat <- function(
	yvls, xvls, dvls = NULL, xaug = NULL, krnl = "ba", band = "and", long = NULL
) {
	lgth <- length(yvls)
	if (is.null(xaug)) {
		xaug <- xvls
	}
	rsds <- qr.resid(qr(cbind(dvls, xaug)), yvls)
	#
	if (is.null(long)) {
		xyvl <- cbind(xvls, yvls)
		wvls <- cbind(xyvl[1:(lgth - 1), ], dvls[2:lgth, ])
		rols <- qr.resid(qr(wvls), xyvl[2:lgth, , drop = FALSE])
		#
		long <- lrvar(rols, krnl, band)$longvar
		lrso <- safesolve(long[-1, -1, drop = FALSE], long[-1, 1, drop = FALSE])
		long <- long[1, 1] - long[1, -1] %*% lrso
	}
	#
	return((lgth * long)/(sum(rsds^2)/lgth))
}