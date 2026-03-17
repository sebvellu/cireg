shinlsstat <- function(
	yvls, zvls, xvls, lale = "aic", krnl = "ba", band = "and",
    dtrn = NULL, step = 2, maxl = NULL, symm = TRUE
) {
    rslt <- dlsfit(
		yvls = yvls, 
		zvls = zvls,
		xvls = xvls,
		lale = lale, 
		rmat = NULL, 
		rvec = NULL,
		wght = NULL,
		dtrn = dtrn,
		step = step,
		maxl = maxl,
		symm = symm
	)
    long <- lrvarstd(yvls, zvls, rslt$vvls, krnl, band)$clrv
	return(drop(kpssstatint(rslt$rsds, long)))
}
