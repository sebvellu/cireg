dcmpr <- function(
	yvls, zpow, xvls, lale = "aic", rmat = NULL, rvec = NULL, wght = NULL,
	krnl = "ba", band = "and", dtrn = NULL, step = 2, maxl = NULL, symm = TRUE
) {
	zvls <- getzvls(zpow, xvls)
	rslt <- dls(
		yvls = yvls,
		zvls = zvls,
		xvls = xvls, 
		lale = lale,
		rmat = rmat,
		rvec = rvec,
		wght = wght,
		krnl = krnl,
		band = band,
		dtrn = dtrn,
		step = step,
		maxl = maxl,
		symm = symm
	)
    rslt$zpow <- zpow
    attr(rslt$zpow, "xnum") <- ncol(xvls)
    rslt$ynum <- NCOL(yvls)
    return(rslt)
}