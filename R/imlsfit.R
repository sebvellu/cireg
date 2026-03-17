
imlsfit <- function(
    yvls, zvls, xvls, zadd = NULL, rmat = NULL, rvec = NULL, wght = NULL
) {
	yvls <- as.matrix(yvls)
	zvls <- as.matrix(zvls)
	xvls <- as.matrix(xvls)
	#
	lgth <- nrow(yvls)
	#ynum <- ncol(yvls)
	znum <- ncol(zvls)
	#xnum <- ncol(xvls)
    #
    s1yv <- apply(yvls, 2, cumsum)
	wvls <- cbind(apply(zvls, 2, cumsum), zadd, xvls)
    #
    rslt <- lsfit(s1yv, wvls, rmat, rvec, wght)
    ufit <- tcrossprod(zvls, rslt$cffs[, seq_len(znum), drop = FALSE])
	ursd <- yvls - ufit
	rfit <- tcrossprod(zvls, rslt$cffr[, seq_len(znum), drop = FALSE])
	rrsd <- yvls - rfit
    #
	s1wv <- apply(wvls, 2, cumsum)
	cvls <- sweep(-s1wv[1:(lgth - 1), , drop = FALSE], 2, s1wv[lgth, ], "+")
	cvls <- rbind(s1wv[lgth, ], cvls)
    #
	rsdm <- qr.resid(qr(cbind(wvls, apply(cvls, 2, cumsum))), s1yv)
    #
    return(list(
        iobs = seq_len(lgth),
        #
        s1yv = s1yv,
        wvls = wvls,
        qrwv = rslt$qrxv,
        #
        cffs = rslt$cffs,
        fits = rslt$fits,
        rsds = rslt$rsds,
        #
        ahlf = rslt$ahlf,
        cfff = rslt$cfff,
        cffr = rslt$cffr,
        fitr = rslt$fitr,
        rsdr = rslt$rsdr,
        #
        ufit = ufit,
        ursd = ursd, 
        rfit = rfit,
        rrsd = rrsd,
        #
        cvls = cvls,
        rsdm = rsdm
    ))
}