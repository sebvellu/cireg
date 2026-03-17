dolslaglead <- function(
	yvls, zvls, xvls, ynum, znum, xnum, maxl, symm = TRUE, lale = "aic",
	dtrn = NULL, step = 2
) {
	info <- Inf
	lvls <- NULL
	#
	if (symm) {
		maxl <- min(maxl)
		for (lags in 0:maxl) {
			rsds <- dlsfit(
				yvls = yvls,
				zvls = zvls,
				xvls = xvls,
				lale = c(lags, lags),
				rmat = NULL,
				rvec = NULL,
				wght = NULL,
				dtrn = dtrn,
				step = step,
                maxl = maxl,
                symm = symm
			)$rsds
			smpl <- nrow(rsds)
			#pars <- (xnum * (2 * lags + 1) + znum) * ynum + ynum^2
			#pars <- (xnum * (2 * lags + 1) + znum + ynum) * ynum
			pars <- (xnum * (2 * lags + 1) + znum + (ynum + 1)/2) * ynum
			if (lale == "aic") {
				temp <- 2 * pars
			} else { #if (lale == "bic") {
				temp <- log(smpl) * pars
			}
			temp <- smpl * ldet(crossprod(rsds)/smpl) + temp
			if (temp < info) {
				lvls <- c(lags, lags)
				info <- temp
			}
		}
		return(lvls)
	} else {
		for (lags in 0:maxl[1]) {
			for (lead in 0:maxl[2]) {
				rsds <- dlsfit(
					yvls = yvls,
					zvls = zvls,
					xvls = xvls,
					lale = c(lags, lead),
				    rmat = NULL,
				    rvec = NULL,
				    wght = NULL,
					dtrn = dtrn,
					step = step,
                    maxl = maxl,
                    symm = symm
				)$rsds
				smpl <- nrow(rsds)
				#pars <- (xnum * (lags + lead + 1) + znum) * ynum + ynum^2
				#pars <- (xnum * (lags + lead + 1) + znum + ynum) * ynum
				pars <- (xnum * (2 * lags + 1) + znum + (ynum + 1)/2) * ynum
				if (lale == "aic") {
					temp <- 2 * pars
				} else { #if (lale == "bic") {
					temp <- log(smpl) * pars
				}
				temp <- smpl * ldet(crossprod(rsds)/smpl) + temp
				if (temp < info) {
					lvls <- c(lags, lead)
					info <- temp
				}
			}
		}
		return(lvls)
	}
}
