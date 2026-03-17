#' Lookup critical values for the 
#' Kwiatkowski, Phillips, Schmidt, Shin (1992) test
#'
#' Uses an internal lookup table of Kwiatkowski, Phillips, Schmidt, Shin (1992)
#' critical values stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupkpss92qntl <- function(prob, dpow, tolr = .Machine$double.eps) {
    tabl <- crit_kpss92 # stored in sysdata.rda
	cndn <- (tabl$dpow == dpow)
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