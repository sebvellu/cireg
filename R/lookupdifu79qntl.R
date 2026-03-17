#' Lookup critical values for Dickey and Fuller (1979)
#'
#' Uses an internal lookup table of Dickey and Fuller (1979) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupdifu79qntl <- function(prob, dpow, stat, tolr = .Machine$double.eps) {
	tabl <- crit_difu79  # stored in sysdata.rda
    cndn <- (tabl$stat == stat) & (tabl$dpow == dpow)
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
