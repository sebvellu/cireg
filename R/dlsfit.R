dlsfit <- function(
	yvls, zvls, xvls, lale = "aic", rmat = NULL, rvec = NULL, wght = NULL,
	dtrn = NULL, step = 2, maxl = NULL, symm = TRUE
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
	if (is.character(lale[[1]])) {
		if (is.null(maxl)) {
			temp <- floor(12 * (lgth/100)^(1/4))
			temp <- min(temp, floor((lgth/(20 * xnum * ynum) - 1)/2))
			temp <- max(1, temp)
            #temp <- max(1, min(2, floor(lgth/(20 * xnum * ynum))))
			maxl <- c(temp, temp)
		} else {
            if (length(maxl) == 1) {
                maxl <- c(maxl, maxl)
            }
        }
		lale <- dolslaglead(
			yvls = yvls,
			zvls = zvls,
			xvls = xvls,
			ynum = ynum,
			znum = znum,
			xnum = xnum,
			maxl = maxl,
			symm = symm,
			lale = lale,
			dtrn = dtrn,
			step = step
		)
		lags <- lale[1]
		lead <- lale[2]
	} else {
		lags <- lale[1]
		lead <- lale[2]
	}
	#
	vvls <- helperkit::diffdtrnd(xvls, dtrn, step)
	#
	lmax <- lgth - lead - 1
	bwrd <- matrix(NA_real_, lmax - lags, xnum * lags)
	fwrd <- matrix(NA_real_, lmax - lags, xnum * lead)
	#
	strt <- 1
	for (indx in seq_len(lags)) {
		fnsh <- indx * xnum
		bwrd[, strt:fnsh] <- vvls[indx:(lmax - lags - 1 + indx), , drop = FALSE]
		strt <- fnsh + 1
	}
	strt <- 1
	for (indx in seq_len(lead)) {
		fnsh <- indx * xnum
		fwrd[, strt:fnsh] <- vvls[(lags + indx + 1):(lmax + indx), , drop = FALSE]
		strt <- fnsh + 1
	}
	#
	indx <- (lags + 1):(lgth - lead)
	vend <- lgth - lags - lead 
	wvls <- cbind(zvls[indx[-1], ], bwrd, vvls[indx[-vend], ], fwrd)
	yadj <- yvls[indx[-1], , drop = FALSE]
	#
    rslt <- lsfit(yadj, wvls, rmat, rvec, wght)
    ufit <- tcrossprod(zvls, rslt$cffs[, seq_len(znum), drop = FALSE])
	ursd <- yvls - ufit
	rfit <- tcrossprod(zvls, rslt$cffr[, seq_len(znum), drop = FALSE])
	rrsd <- yvls - rfit
	#
    return(list(
        lale = lale,
        iobs = indx[-1],
        #
        yadj = yadj,
        wvls = wvls,
        vvls = vvls,
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
		rrsd = rrsd
    ))
}