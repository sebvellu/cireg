#' Lookup critical values for Phillips and Ouliaris (1990)
#'
#' Uses an internal lookup table of Phillips and Ouliaris (1990) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupphou90qntl <- function(
	prob, dpow, stat, xnum, tolr = .Machine$double.eps
) {
	tabl <- crit_phou90  # stored in sysdata.rda
    cndn <- (tabl$stat == stat) & (tabl$xnum == xnum) & (tabl$dpow == dpow)
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
