#' Lookup critical values for CT and PU cpr
#'
#' Uses an internal lookup table of CT and PU cpr critical values 
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupwa23qntl <- function(
	prob, dpow, xpow, xnum, stat, tolr = .Machine$double.eps
) {
    tabl <- crit_wa23 # stored in sysdata.rda
	cndn <- (tabl$stat == stat) & (tabl$xnum == xnum) & (tabl$xpow == xpow)
	cndn <- cndn & (tabl$dpow == dpow)
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