fmlsfit <- function(
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
	vvls <- helperkit::diffdtrnd(xvls, dtrn, step)
	lrvr <- lrvarstd(yvls, zvls, vvls, krnl, band)
    lrso <- lrvr$lrso
    clrv <- lrvr$clrv
    lrvr <- lrvr$lrvr
    #
	xidx <- ynum + 1:xnum
	yidx <- 1:ynum
    #
	bico <- lrvr$deltvar[xidx, yidx] - lrvr$deltvar[xidx, xidx] %*% lrso
    bico <- rowSums(dzdx, dims = 2) %*% bico
    #
    ymod <- yvls[-1, , drop = FALSE]
	yadj <- ymod - vvls %*% lrso
	wvls <- zvls[-1, , drop = FALSE]
	#
    qrwv <- qr(wvls)
	rqrw <- qr.R(qrwv)
	#whiv <- helperkit::qr_safe_solve(t(rqrw), diag(1, ncol(rqrw)))
	whiv <- helperkit::safe_forwardsolve(t(rqrw), diag(1, ncol(rqrw)))
	cent <- crossprod(whiv)
    #
    # lsfit
    yols <- yadj - (wvls %*% cent) %*% bico
    #cffs <- qr.coef(qrwv, yadj) - cent %*% bico # ==
    cffs <- qr.coef(qrwv, yols)
    rsds <- yadj - wvls %*% cffs
    fits <- ymod - rsds
    cffs <- t(cffs)
    #
    if (is.null(rmat) || is.null(rvec) || is.null(wght)) {
        ahlf <- NULL
        cfff <- c(t(cffs))
        cffr <- cffs
        fitr <- fits
        rsdr <- rsds
    } else {
	    rslt <- helperkit::rfmls(ymod, yadj, yols, xvls, rmat, rvec, wght)
	    ahlf <- rslt$ahlf
        cfff <- rslt$cfff
	    cffr <- rslt$cffr
	    fitr <- rslt$fitr
	    rsdr <- rslt$rsdr
    }
    #
    ufit <- tcrossprod(zvls, cffs)
    ursd <- yvls - ufit
    rfit <- tcrossprod(zvls, cffr)
    rrsd <- yvls - rfit
	#
    return(list(
        iobs = seq_len(lgth - 1) + 1,
        long = lrvr$longvar,
        clrv = clrv,
        lrso = lrso,
        #
        vvls = vvls,
        #
        ymod = ymod,
        yadj = yadj,
        yols = yols,
        wvls = wvls,
        qrwv = qrwv,
        #
        cffs = cffs,
        fits = fits,
        rsds = rsds,
        #
        ahlf = ahlf,
        cfff = cfff,
        cffr = cffr,
        fitr = fitr,
        rsdr = rsdr,
        #
        ufit = ufit,
        ursd = ursd,
        rfit = rfit,
        rrsd = rrsd
    ))
}