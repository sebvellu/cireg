fmcmpr <- function(
	yvls, zpow, xvls, krnl = "ba", band = "and",
    rmat = NULL, rvec = NULL, wght = NULL, dtrn = NULL, step = 2
) {
	zvls <- getzvls(zpow, xvls)
	dzdx <- getdzdx(zpow, xvls)
	rslt <- fmls(
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
    )
    rslt$zpow <- zpow
    attr(rslt$zpow, "xnum") <- ncol(xvls)
    rslt$ynum <- NCOL(yvls)
    return(rslt)
}
