shincmprstat <- function(
	yvls, zpow, xvls, lale = "aic", krnl = "ba", band = "and",
    dtrn = NULL, step = 2, maxl = NULL, symm = TRUE
) {
    zvls <- getzvls(zpow, xvls)
	return(shinlsstat(
        yvls = yvls, 
        zvls = zvls,
        xvls = xvls,
        lale = lale,
        krnl = krnl,
        band = band,
        dtrn = dtrn,
        step = step, 
        maxl = maxl,
        symm = symm
    ))
}
