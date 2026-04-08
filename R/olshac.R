#' @export

olshac <- function(
	yvls, dvls, xvls, krnl = "ba", band = "and",
    rmat = NULL, rvec = NULL, wght = NULL, dtrn = NULL, step = 2
) {
	zvls <- cbind(dvls, xvls)
	return(ls(
        yvls = yvls,
        zvls = zvls,
        xvls = xvls, 
        krnl = krnl,
        band = band,
        rmat = rmat,
        rvec = rvec,
        wght = wght, 
        dtrn = dtrn,
        step = step
    ))
}
