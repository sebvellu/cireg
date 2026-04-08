
lrvarstd <- function(yvls, zvls, vvls, krnl = "ba", band = "band") {
	yvls <- as.matrix(yvls)
	#
	rols <- cbind(qr.resid(qr(zvls), yvls)[-1, ], vvls)
	#rols <- yvls - zvls %*% qr_safe_solve(zvls, yvls)
	#rols <- cbind(rols[-1, ], vvls)
    lrvr <- lrvar::lrvar(rols, krnl, band)
	long <- lrvr$longvar
	#
	xidx <- ncol(yvls) + seq_len(ncol(vvls))
	yidx <- seq_len(ncol(yvls))
	#
	lrso <- helperkit::qr_safe_solve(
		mata = long[xidx, xidx, drop = FALSE], 
		matb = long[xidx, yidx, drop = FALSE]
	)
	clrv <- long[yidx, yidx] - long[yidx, xidx] %*% lrso
    #
    return(list(
        rols = rols,
        lrvr = lrvr,
        lrso = lrso,
        clrv = clrv
    ))
}