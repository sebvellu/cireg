rfmls <- function(yvls, yadj, yols, xvls, rmat, rvec, wght) {
	whlf <- safe_chol(wght)
	ahlf <- (whlf %x% xvls) %*% rmat
	yadj <- c(tcrossprod(yols, whlf)) - (whlf %x% xvls) %*% rvec
	#
	cfff <- qr_safe_solve(ahlf, yadj)
	cffr <- rmat %*% cfff + rvec
	cffr <- matrix(cffr, ncol(xvls), NCOL(yvls))
	#
	rsdr <- yadj - xvls %*% cffr
    fitr <- yvls - rsdr
	cffr <- t(cffr)
    #
    return(list(
        ahlf = ahlf,
        cfff = cfff,
        cffr = cffr,
        fitr = fitr,
        rsdr = rsdr
    ))
}
