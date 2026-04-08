imcmprfit <- function(
	yvls, zpow, xvls, zadd = NULL, rmat = NULL, rvec = NULL, wght = NULL
) {
	zvls <- getzvls(zpow, xvls)
	rslt <- imlsfit(
		yvls = yvls,
		zvls = zvls,
		xvls = xvls,
		zadd = zadd,
		rmat = rmat,
		rvec = rvec,
		wght = wght
	)
    rslt$zpow <- zpow
    attr(rslt$zpow, "xnum") <- ncol(xvls)
    rslt$ynum <- NCOL(yvls)
	return(rslt)
}