imls <- function(
	yvls, zvls, xvls, zadd = NULL, rmat = NULL, rvec = NULL, wght = NULL, 
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
	rslt <- imlsfit(yvls, zvls, xvls, zadd, NULL, NULL, NULL)
    s1yv <- rslt$s1yv
    wvls <- rslt$wvls
    qrwv <- rslt$qrwv
	cffs <- rslt$cffs
	fits <- rslt$fits
	rsds <- rslt$rsds
	ufit <- rslt$ufit
	ursd <- rslt$ursd
    #
	cvls <- rslt$cvls
	rsdm <- rslt$rsdm
	#
	wnum <- ncol(wvls)
	#
	# rqrw <- qr.R(qrwv)
	# #whiv <- qr_safe_solve(rqrw, qr_safe_solve(t(rqrw), diag(1, ncol(rqrw))))
	# temp <- safe_forwardsolve(t(rqrw), diag(1, ncol(rqrw)))
	# whiv <- safe_backsolve(rqrw, temp)
	# cent <- cvls %*% whiv
	# cent <- crossprod(cent)
	cent <- crossprod(cvls %*% sqrtinvprod(qrwv, TRUE))
	#
	# Standard asymptotic theory, covariance matrix of vec(cffs')
	#
	lrvr <- lrvarstd(yvls, zvls, diffdtrnd(xvls, dtrn, step), krnl, band)
	long <- lrvr$longvar
	clrv <- lrvr$clrv
	#
	cvrs <- clrv %x% cent
	csds <- matrix(sqrt(diag(cvrs)), wnum, ynum)
	csds <- t(csds)
	#
	# Fixed-b asymptotic theory:
	#
	lrvr <- lrvar(rbind(rsdm[1, , drop = FALSE], diff(rsdm)), krnl, band)
	clrb <- lrvr$longvar
	fixb <- lrvr$bwdh/lgth
	#
	cvrb <- clrb %x% cent
	csdb <- matrix(sqrt(diag(cvrb)), wnum, ynum)
	csdb <- t(csdb)
	#
	# Restricted Case:
	#
	if (is.null(rmat)) {
		rmat <- diag(1, (znum + xnum) * ynum)
	}
	if (is.null(rvec)) {
		rvec <- numeric((znum + xnum) * ynum)
	}
	if (is.null(wght)) {
		wght <- diag(1, ynum)
	} else if (identical(wght, "PaOg1")) {
		wght <- safe_inv(long[1:ynum, 1:ynum])
	} else if (identical(wght, "PaOg2")) {
		wght <- safe_inv(clrv)
	}
	#
	rslt <- rls(s1yv, wvls, rmat, rvec, wght)
	ahlf <- rslt$ahlf
	cfff <- rslt$cfff
	cffr <- rslt$cffr
	fitr <- rslt$fitr
	rsdr <- rslt$rsdr
	#
	rfit <- tcrossprod(zvls, cffr[, seq_len(znum), drop = FALSE])
	rrsd <- yvls - rfit
	#
	# rqra <- qr.R(qr(ahlf))
	# #ahiv <- qr_safe_solve(rqra, qr_safe_solve(t(rqra), diag(1, ncol(rqra))))
	# temp <- safe_forwardsolve(t(rqra), diag(1, ncol(rqra)))
	# ahiv <- safe_backsolve(rqra, temp)
	ahiv <- sqrtinvprod(qr(ahlf), TRUE)
	#
	#bhlf <- ((safe_chol(clrv) %*% wght) %x% cvls) %*% rmat
	bhlf <- calcbhlf(cvls, clrv, rmat, wght)
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
		rsdm = rsdm,
		fixb = fixb,
		cvrb = cvrb,
		csdb = csdb,
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
