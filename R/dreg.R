dreg <- function(
	yvls, dvls, xvls, lale = "aic", rmat = NULL, rvec = NULL, wght = NULL,
	krnl = "ba", band = "and", dtrn = NULL, step = 2, maxl = NULL, symm = TRUE
) {
	zvls <- cbind(dvls, xvls)
	return(dls(
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
	))
}