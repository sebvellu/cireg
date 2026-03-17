rls <- function(yvls, xvls, rmat, rvec, wght) {
	whlf <- safe_chol(wght)
	ahlf <- (whlf %x% xvls) %*% rmat
	yadj <- c(tcrossprod(yvls, whlf)) - (whlf %x% xvls) %*% rvec
	#
	cfff <- qr_safe_solve(ahlf, yadj)
	cffr <- rmat %*% cfff + rvec
	cffr <- matrix(cffr, ncol(xvls), NCOL(yvls))
	#
	fitr <- xvls %*% cffr
	rsdr <- yvls - fitr
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
