getpovrcmprqntl <- function(
    prob, zpow, lgth = 1000, simu = 10000, tolr = .Machine$double.eps
) {
	if (!is.matrix(zpow)) {
		if (length(zpow) == 1) {
			temp <- try(lookupphou90qntl(
                prob = prob,
                dpow = zpow,
                stat = "pu",
                xnum = attr(zpow, "xnum"),
                tolr = tolr
            ), TRUE)
			if (!inherits(temp, "try-error")) {
				if (!any(is.na(temp))) {
					return(temp)
				}
			}
			zpow <- c(zpow, rep(1, attr(zpow, "xnum")))
		}
		dpow <- zpow[1]
		xpow <- zpow[-1]
    	if (sum(xpow == 1) > length(xpow) - 2) {
            temp <- try(lookupwa23qntl(
                prob = prob,
                dpow = dpow,
                xpow = max(xpow),
                xnum = length(xpow),
                stat = "pu",
                tolr = tolr
            ), TRUE)
            if (!inherits(temp, "try-error")) {
                if (!any(is.na(temp))) {
                    return(temp)
                }
            }
	    }
		zpow <- cprtocmprpow(dpow, xpow)
	}
	return(povrcmprquant(prob, zpow, lgth, simu))
}