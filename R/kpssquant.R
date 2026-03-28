kpssquant <- function(
	prob, dpow = NULL, lgth = 1000, simu = 10000
) {
    rslt <- numeric(simu)
	#
	for (indx in 1:simu) {
        tsrs <- stats::rnorm(lgth)
	    if (is.null(dpow)) {
            rsds <- tsrs
	    } else {
            dvls <- polydet(lgth, dpow)
		    rsds <- qr.resid(qr(dvls), tsrs)
	    }
        rslt[indx] <- kpssstatint(rsds, 1)
    }
    return(get_quantile(rslt, prob))
}
