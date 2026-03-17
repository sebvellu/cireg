#' Lookup critical values for Vogeland and Wagner (2013)
#'
#' Uses an internal lookup table of Vogeland and Wagner (2013) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupvowa13qntl <- function(
	prob, dpow, step, stat, krnl, bfrc, tolr = .Machine$double.eps
) {
    tabl <- crit_vowa13 # stored in sysdata.rda
	cndn <- (tabl$stat == stat) & (tabl$krnl == krnl) & (tabl$step == step)
	cndn <- cndn & (tabl$dpow == dpow) & (abs(tabl$bfrc - bfrc) < tolr)
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
