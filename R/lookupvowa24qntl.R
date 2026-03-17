#' Lookup critical values for Vogelsang and Wagner (2024)
#'
#' Uses an internal lookup table of Vogelsang and Wagner (2024) critical values
#' stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupvowa24qntl <- function(
	prob, zpow, hypo, krnl, bfrc, tolr = .Machine$double.eps
) {
    tabl <- crit_vowa24  # stored in sysdata.rda
    znum <- ncol(zpow)
    zpow <- paste(sort(apply(zpow, 1, cantor_tuple)), collapse = "-")
	cndn <- (tabl$krnl == krnl) & (tabl$znum == znum) & (tabl$zpow == zpow)
	cndn <- cndn & (tabl$hypo == hypo) & (abs(tabl$bfrc - bfrc) < tolr)
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
