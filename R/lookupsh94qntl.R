#' Lookup critical values for the Shin (1994) test
#'
#' Uses an internal lookup table of Shin (1994) critical values stored 
#' in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupsh94qntl <- function(prob, xnum, dpow, tolr = .Machine$double.eps) {
    tabl <- crit_sh94 # stored in sysdata.rda
	cndn <- (tabl$xnum == xnum) & (tabl$dpow == dpow)
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