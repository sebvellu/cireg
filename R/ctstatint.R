ctstatint <- function(s1rs, long = NULL) {
    lgth <- nrow(s1rs)
    stat <- numeric(3)
    #
	vala <- crossprod(s1rs)
	valb <- crossprod(apply(s1rs, 2, cumsum))
    if (is.null(long)) {
	    stat[1] <- sum(diag(vala))/lgth^2
	    stat[2] <- sum(diag(valb))/lgth^4
    } else {
	    stat[1] <- sum(diag(helperkit::safesolve(long, vala)))/lgth^2
	    stat[2] <- sum(diag(helperkit::safesolve(long, valb)))/lgth^4	
    }
    stat[3] <- sum(diag(helperkit::safesolve(vala, valb)))/lgth^2
    #
	return(stat)
}