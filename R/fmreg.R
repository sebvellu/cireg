#' @export
#' 
fmreg <- function(
	yvls, dvls, xvls, krnl = "ba", band = "and",
    rmat = NULL, rvec = NULL, wght = NULL, dtrn = NULL, step = 2
) {
	zvls <- cbind(dvls, xvls)
	dzdx <- rbind(matrix(0, ncol(dvls), ncol(xvls)), diag(ncol(xvls)))
	dzdx <- replicate(NROW(yvls), dzdx)
	return(fmls(
        yvls = yvls,
        zvls = zvls,
        xvls = xvls, 
        dzdx = dzdx,
        krnl = krnl,
        band = band,
        rmat = rmat,
        rvec = rvec,
        wght = wght, 
        dtrn = dtrn,
        step = step
    ))
}