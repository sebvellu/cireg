dregfit <- function(
	yvls, dvls, xvls, lale = "aic", rmat = NULL, rvec = NULL, wght = NULL,
	dtrn = NULL, step = 2, maxl = NULL, symm = TRUE
) {
	zvls <- cbind(dvls, xvls)
	return(dlsfit(
		yvls = yvls,
		zvls = zvls,
		xvls = xvls, 
		lale = lale,
		rmat = rmat,
		rvec = rvec,
		wght = wght,
		dtrn = dtrn,
		step = step,
		maxl = maxl,
		symm = symm
	))
}