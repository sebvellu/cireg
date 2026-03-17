#' Lookup critical values for Veldhuis and Wagner (????) CT statitics
#'
#' Uses an internal lookup table of Veldhuis and Wagner (????) CT 
#' critical values stored in `sysdata.rda`.
#' 
#' @keywords internal
#' 
lookupvewa24ctqntl <- function(
    prob, ynum, zpow, stat, tolr = .Machine$double.eps
) {
    tabl <- crit_vewa24_ct # stored in sysdata.rda
    znum <- ncol(zpow)
    zpow <- paste(sort(apply(zpow, 1, cantor_tuple)), collapse = "-")
	cndn <- (tabl$ynum == ynum) & (tabl$zpow == zpow) & (tabl$znum == znum)
    cndn <- cndn & (tabl$stat == stat)
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
