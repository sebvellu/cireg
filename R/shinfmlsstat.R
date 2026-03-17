shinfmlsstat <- function(
	yvls, zvls, xvls, dzdx, krnl = "ba", band = "and", dtrn = NULL, step = 2
) {
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
	return(drop(kpssstatint(rslt$rsds, rslt$clrv)))
}
