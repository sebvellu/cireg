shinfmcmprstat <- function(
	yvls, zpow, xvls, krnl = "ba", band = "and", dtrn = NULL, step = 2
) {
	zvls <- getzvls(zpow, xvls)
	dzdx <- getdzdx(zpow, xvls)
	return(shinfmlsstat(
        yvls = yvls, 
        zvls = zvls,
        xvls = xvls,
        dzdx = dzdx,
        krnl = krnl,
        band = band,
        dtrn = dtrn,
        step = step
    ))
}
