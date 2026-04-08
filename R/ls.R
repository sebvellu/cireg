ls <- function(
	yvls, zvls, xvls, rmat = NULL, rvec = NULL, wght = NULL,
	krnl = "ba", band = "and", dtrn = NULL, step = 2
) {
	yvls <- as.matrix(yvls)
	zvls <- as.matrix(zvls)
	xvls <- as.matrix(xvls)
	#
	lgth <- nrow(yvls)
	ynum <- ncol(yvls)
	znum <- ncol(zvls)
	xnum <- ncol(xvls)
	#
    rslt <- lsfit(
		yvls = yvls, 
		xvls = zvls,
		rmat = NULL, 
		rvec = NULL,
		wght = NULL
	)
    wvls <- rslt$xvls
    qrwv <- rslt$qrxv
	cffs <- rslt$cffs
	fits <- rslt$fits
	rsds <- rslt$rsds
	#ufit <- rslt$ufit
	#ursd <- rslt$ursd
    #
    cvls <- wvls
	wnum <- ncol(wvls)
	#
	cent <- crossprod(helperkit::sqrtinvprod(qrwv, FALSE))
	#
	# Standard asymptotic theory, covariance matrix of vec(cffs')
	#
	vvls <- helperkit::diffdtrnd(xvls, dtrn, step)
	lrvr <- lrvarstd(yvls, zvls, vvls, krnl, band)
	long <- lrvr$longvar
	clrv <- lrvr$clrv
	#
	cvrs <- clrv %x% cent
	csds <- matrix(sqrt(diag(cvrs)), wnum, ynum)
	csds <- t(csds)
	#
	# Restricted Case:
	#
	if (is.null(rmat)) {
		rmat <- diag(1, wnum * ynum)
	}
	if (is.null(rvec)) {
		rvec <- numeric(wnum * ynum)
	}
	if (is.null(wght)) {
		wght <- diag(1, ynum)
	} else if (identical(wght, "PaOg1")) {
		wght <- helperkit::safe_inv(long[1:ynum, 1:ynum])
	} else if (identical(wght, "PaOg2")) {
		wght <- helperkit::safe_inv(clrv)
	}
	#
	rslt <- helperkit::rls(yvls, wvls, rmat, rvec, wght)
	ahlf <- rslt$ahlf
	cfff <- rslt$cfff
	cffr <- rslt$cffr
	fitr <- rslt$fitr
	rsdr <- rslt$rsdr
	#
	#rfit <- tcrossprod(zvls, cffr)
	#rrsd <- yvls - rfit
	#
    ahiv <- helperkit::sqrtinvprod(qr(ahlf), TRUE)
    #
    #bhlf <- ((helperkit::safe_chol(clrv) %*% wght) %x% cvls) %*% rmat # cvls == wvls
	bhlf <- helperkit::calcbhlf(cvls, clrv, rmat, wght) # cvls == wvls
    #
	cvrf <- bhlf %*% ahiv
	cvrr <- tcrossprod(cvrf, rmat)
	cvrr <- crossprod(cvrr)
	cvrf <- crossprod(cvrf)
	#
	csdr <- matrix(sqrt(diag(cvrr)), wnum, ynum)
	csdr <- t(csdr)
	#
	return(list(
		cffs = cffs, #cffu
		fits = fits, #fitu
		rsds = rsds, #rsdu
		cvrs = cvrs, #cvru
		csds = csds, #csdu
		clrv = clrv,
		cffr = cffr,
		fitr = fitr,
		rsdr = rsdr,
		cvrr = cvrr, 
		csdr = csdr,
		cfff = cfff,
		cvrf = cvrf,
		ufit = fits,
		ursd = rsds,
		rfit = fitr,
		rrsd = rsdr
	))
}
