ctimlsstat <- function(
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
	vvls <- helperkit::diffdtrnd(xvls, dtrn, step)
    long <- lrvarstd(yvls, zvls, vvls, krnl, band)$clrv
    s1rs <- rslt$rsds
	return(ctstatint(s1rs, long))
}