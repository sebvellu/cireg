imreg <- function(
	yvls, dvls, xvls, zadd = NULL, rmat = NULL, rvec = NULL, wght = NULL, 
	krnl = "ba", band = "and", dtrn = NULL, step = 2
) {
	zvls <- cbind(dvls, xvls)
	return(imls(
		yvls = yvls,
		zvls = zvls,
		xvls = xvls,
		zadd = zadd,
		rmat = rmat,
		rvec = rvec,
		wght = wght,
		krnl = krnl,
		band = band,
		dtrn = dtrn,
		step = step
	))
}