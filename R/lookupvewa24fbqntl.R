#' Lookup critical values for Veldhuis and Wanger (...)
#'
#' Uses an internal lookup table of Veldhuis and Wagner (...) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupvewa24fbqntl <- function(
	prob, ynum, zpow, stat, hypo, krnl, bfrc, tolr = .Machine$double.eps
) {
    tabl <- crit_vewa24_fb   # stored in sysdata.rda
    znum <- ncol(zpow)
    zpow <- paste(sort(apply(zpow, 1, cantor_tuple)), collapse = "-")
	cndn <- (tabl$stat == stat) & (tabl$krnl == krnl) & (tabl$ynum == ynum)
	cndn <- cndn & (tabl$znum == znum) & (tabl$zpow == zpow)
    cndn <- cndn & (tabl$hypo == hypo)
	cndn <- cndn & (abs(tabl$bfrc - bfrc) < tolr)
	rslt <- numeric(length(prob))
	for (pidx in seq_along(prob)) {
		temp <- cndn & (abs(tabl$prob - prob[pidx]) < tolr)
		if (any(temp)) {
			rslt[pidx] <- tabl[temp, 2]
		} else {
			rslt[pidx] <- NA_real_
		}
	}
	#names(rslt) <- prob
	return(rslt)
}
