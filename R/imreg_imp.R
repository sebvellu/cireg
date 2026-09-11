imreg_imp <- function(
	yvls, dvls, xvls, clrv, long, zadd = NULL, rmat = NULL, rvec = NULL,
	wght = NULL
) {
	zvls <- cbind(dvls, xvls)
	return(imls_imp(
		yvls = yvls,
		zvls = zvls,
		xvls = xvls,
		clrv = clrv,
		long = long,
		zadd = zadd,
		rmat = rmat,
		rvec = rvec,
		wght = wght
	))
}