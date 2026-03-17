#' Lookup critical values for Veldhuis and Wagner (????)
#'
#' Uses an internal lookup table of Veldhuis and Wagner (????) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupvewa23qntl <- function(
	prob, dpow, xnum, stat, krnl, bfrc, tolr = .Machine$double.eps
) {
    tabl <- crit_vewa23 # stored in sysdata.rda
	cndn <- (tabl$stat == stat) & (tabl$krnl == krnl) & (tabl$xnum == xnum)
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
