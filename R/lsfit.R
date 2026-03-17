lsfit <- function(yvls, xvls, rmat = NULL, rvec = NULL, wght = NULL) {
	lgth <- nrow(yvls)
	qrxv <- qr(xvls)
    #
	#cffs <- qr_safe_solve(xvls, yvls)
	#fits <- xvls %*% cffs
	#rsds <- yvls - fits
	#cffs <- t(cffs)
    cffs <- qr.coef(qrxv, yvls)
    fits <- qr.fitted(qrxv, yvls)
    rsds <- qr.resid(qrxv, yvls) #yvls - fits
    cffs <- t(cffs)
    #
    if (is.null(rmat) || is.null(rvec) || is.null(wght)) {
        ahlf <- NULL
        cfff <- c(t(cffs))
        cffr <- cffs
        fitr <- fits
        rsdr <- rsds
    } else {
	    rslt <- rls(yvls, xvls, rmat, rvec, wght)
	    ahlf <- rslt$ahlf
        cfff <- rslt$cfff
	    cffr <- rslt$cffr
	    fitr <- rslt$fitr
	    rsdr <- rslt$rsdr
    }
    #
    return(list(
        yvls = yvls,
        xvls = xvls,
        qrxv = qrxv,
        #
        cffs = cffs,
        fits = fits,
        rsds = rsds,
        #
        ahlf = ahlf,
        cfff = cfff,
        cffr = cffr,
        fitr = fitr,
        rsdr = rsdr
    ))
}