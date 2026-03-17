computelrmat <- function(tsrs, wght) { # tsrs = tsrs/sqrt(lgth)
	lgth <- nrow(tsrs)
	nums <- ncol(tsrs)
	#
	half <- matrix(0, nums, nums) #half <- 0
	for (h in seq_along(wght)) {
		temp <- crossprod(
			tsrs[1:(lgth - h), , drop = FALSE],
			tsrs[(h + 1):lgth, , drop = FALSE]
		)
		half <- half + wght[h]  * temp
	}
	shrt <- crossprod(tsrs)
	delt <- shrt + half
	long <- delt + t(half)
	return(list(shrt = shrt, half = half, delt = delt, long = long))
}
