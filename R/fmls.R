fmls <- function(
	yvls, zvls, xvls, dzdx, krnl = "ba", band = "and",
    rmat = NULL, rvec = NULL, wght = NULL, dtrn = NULL, step = 2
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
    rslt <- fmlsfit(
		yvls = yvls, 
		zvls = zvls,
		xvls = xvls,
		dzdx = dzdx,
        krnl = krnl,
        band = band,
		rmat = NULL, 
		rvec = NULL,
		wght = NULL,
		dtrn = dtrn,
		step = step
	)
	vvls <- rslt$vvls
	long <- rslt$long
	clrv <- rslt$clrv
    ymod <- rslt$ymod
	yadj <- rslt$yadj
    yols <- rslt$yols
	wvls <- rslt$wvls
	qrwv <- rslt$qrwv
	cffs <- rslt$cffs
	fits <- rslt$fits
	rsds <- rslt$rsds
	ufit <- rslt$ufit
	ursd <- rslt$ursd
	#
	cvls <- wvls
	wnum <- ncol(wvls)
	#
	#rqrw <- qr.R(qrwv)
	##whiv <- qr_safe_solve(rqrw, qr_safe_solve(t(rqrw), diag(1, ncol(rqrw))))
	#temp <- safe_forwardsolve(t(rqrw), diag(1, ncol(rqrw)))
	#whiv <- safe_backsolve(rqrw, temp)
	#cent <- cvls %*% whiv
	#cent <- crossprod(cent)
	#
	# rqrw <- qr.R(qrwv)
	# #whiv <- qr_safe_solve(t(rqrw), diag(1, ncol(rqrw)))
	# whiv <- safe_forwardsolve(t(rqrw), diag(1, ncol(rqrw)))
	# cent <- crossprod(whiv)
    cent <- crossprod(helperkit::sqrtinvprod(qrwv, FALSE)) # ==
    #cent <- crossprod(cvls %*% helperkit::sqrtinvprod(qrwv, TRUE))
	#
	#lrvr <- lrvarstd(ymod, zvls, diffdtrnd(xvls, dtrn, step), krnl, band)
	#long <- lrvr$longvar
	#clrv <- lrvr$clrv
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
	rslt <- helperkit::rfmls(ymod, yadj, yols, wvls, rmat, rvec, wght)
	ahlf <- rslt$ahlf
	cfff <- rslt$cfff
	cffr <- rslt$cffr
	fitr <- rslt$fitr
	rsdr <- rslt$rsdr
	#
	rfit <- tcrossprod(zvls, cffr)
	rrsd <- yvls - rfit
	#
	#rqra <- qr.R(qr(ahlf))
	##ahiv <- helperkit::qr_safe_solve(rqra, qr_safe_solve(t(rqra), diag(1, ncol(rqra))))
	#temp <- helperkit::safe_forwardsolve(t(rqra), diag(1, ncol(rqra)))
	#ahiv <- helperkit::safe_backsolve(rqra, temp)
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
		ufit = ufit,
		ursd = ursd,
		rfit = rfit,
		rrsd = rrsd
	))
}