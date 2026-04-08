imregfit <- function(
	yvls, dvls, xvls, zadd = NULL, rmat = NULL, rvec = NULL, wght = NULL
) {
	zvls <- cbind(dvls, xvls)
	return(imlsfit(
		yvls = yvls,
		zvls = zvls,
		xvls = xvls,
		zadd = zadd,
		rmat = rmat,
		rvec = rvec,
		wght = wght
	))
}