getkpssqntl <- function(
    prob, dpow, lgth = 1000, simu = 10000, tolr = .Machine$double.eps
) {
	if (!is.null(dpow)) {
		temp <- try(lookupkpss92qntl(
			prob = prob,
			dpow = dpow,
			tolr = tolr
		), TRUE)
		if (!inherits(temp, "try-error")) {
			if (!any(is.na(temp))) {
				return(temp)
			}
		}
	}
	return(kpssquant(prob, dpow, lgth, simu))
}