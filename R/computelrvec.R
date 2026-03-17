computelrvec <- function(tsrs, wght) { # tsrs = tsrs/sqrt(lgth)
	lgth <- length(tsrs)
	#
	half <- 0
	for (h in seq_along(wght)) {
		temp <- sum(tsrs[1:(lgth - h)] * tsrs[(h + 1):lgth])
		half <- half + wght[h] * temp
	}
	shrt <- sum(tsrs^2)
	delt <- shrt + half
	long <- delt + half
	return(list(shrt = shrt, half = half, delt = delt, long = long))
}
