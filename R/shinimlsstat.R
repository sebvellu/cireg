shinimlsstat <- function(
	yvls, zvls, xvls, zadd = NULL, krnl = "ba", band = "and",
    dtrn = NULL, step = 2
) {
    rslt <- imlsfit(
		yvls = yvls, 
		zvls = zvls,
		xvls = xvls,
        zadd = zadd,
		rmat = NULL, 
		rvec = NULL,
		wght = NULL
	)
    long <- lrvarstd(yvls, zvls, diffdtrnd(xvls, dtrn, step), krnl, band)$clrv
    rsds <- c(rslt$rsds[1], diff(rslt$rsds))
	return(drop(kpssstatint(rsds, long)))
}
