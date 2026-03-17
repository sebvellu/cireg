difustat <- function(yvls, dvls = NULL, xvls = NULL, step = 2) {
	lgth <- length(yvls)
	if (!is.null(xvls)) {
		yvls <- qr.resid(qr(cbind(dvls, xvls)), yvls)
		dvls <- NULL
	}
	if (is.null(dvls)) {
		ylag <- yvls[-lgth]
		yvls <- yvls[-1]
	} else {
		if (step == 1) {
			qrdv <- qr(dvls[-1, , drop = FALSE])
			ylag <- qr.resid(qrdv, yvls[-lgth])
			yvls <- qr.resid(qrdv, yvls[-1])
		} else { # if (step == 2) {
			yvls <- qr.resid(qr(dvls), yvls)
			ylag <- yvls[-lgth]
			yvls <- yvls[-1]
		}
	}
	#
	deno <- sum(ylag^2)
	coef <- sum(yvls * ylag)/deno
	rsds <- yvls - coef * ylag
	#
	shrt <- sum(rsds^2)/lgth
    cffs <- lgth * (coef - 1)
	tsta <- (coef - 1)/sqrt(shrt/deno)
	#
	return(list(
        lgth = lgth,
        yvls = yvls,
        ylag = ylag,
        deno = deno,
		coef = coef,
		rsds = rsds,
		shrt = shrt,
        cffs = cffs,
		tsta = tsta
	))
}
