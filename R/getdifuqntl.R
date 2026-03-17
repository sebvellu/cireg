getdifuqntl <- function(
	prob, dpow, xnum = 0, lgth = 1000, simu = 10000, tolr = .Machine$double.eps
) {
	if (xnum == 0) {
		tema <- try(lookupdifu79qntl(
			prob = prob,
			dpow = dpow,
			stat = "zcff",
			tolr = tolr
		), TRUE)
	} else {
		tema <- try(lookupphou90qntl(
			prob = prob,
			dpow = dpow,
			stat = "zcff",
			xnum = xnum,
			tolr = tolr
		), TRUE)
	}
	if (!inherits(tema, "try-error")) {
		if (!any(is.na(tema))) {
			if (xnum == 0) {
				temb <- lookupdifu79qntl(
					prob = prob,
					dpow = dpow,
					stat = "ztst",
					tolr = tolr
				)
			} else {
				temb <- lookupphou90qntl(
					prob = prob,
					dpow = dpow,
					stat = "ztst",
					xnum = xnum,
					tolr = tolr
				)
			}
			return(drop(cbind(tema, temb, deparse.level = 0)))
		}
	}
	return(difuquant(prob, dpow, xnum, 1, lgth, simu))
}