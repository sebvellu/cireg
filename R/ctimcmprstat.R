ctimcmprstat <- function(
	yvls, zpow, xvls, zadd = NULL, krnl = "ba", band = "and",
    dtrn = NULL, step = 2
) {
	zvls <- getzvls(zpow, xvls)
	return(ctimlsstat(
        yvls = yvls, 
        zvls = zvls,
        xvls = xvls,
        zadd = zadd,
        krnl = krnl,
        band = band,
        dtrn = dtrn,
        step = step
    ))
}
