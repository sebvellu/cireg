#' Lookup critical values for Vogelsang and Wagner (2014)
#'
#' Uses an internal lookup table of Vogelsang and Wagner (2014) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupvowa14qntl <- function(
	prob, dpow, stat, xnum, hypo, krnl, bfrc, tolr = .Machine$double.eps
) {
    tabl <- crit_vowa14   # stored in sysdata.rda
	cndn <- (tabl$stat == stat) & (tabl$krnl == krnl) & (tabl$hypo == hypo)
	cndn <- cndn & (tabl$xnum == xnum) & (tabl$dpow == dpow)
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