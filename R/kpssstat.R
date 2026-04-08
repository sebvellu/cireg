kpssstat <- function(tsrs, dvls = NULL, krnl = "ba", band = "and") { # no fixed-b variant #match.fun
	if (is.null(dvls)) {
        rsds <- tsrs
	} else {
		rsds <- qr.resid(qr(dvls), tsrs)
	}
	long <- lrvar::lrvar(rsds, krnl, band)$longvar
    return(kpssstatint(rsds, long))
}