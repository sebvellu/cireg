getshinimcmprqntl <- function(
	prob, zpow, zadd, lgth = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    if (is.null(zadd)) {
        if (!is.matrix(zpow)) {
            if (length(zpow) == 1) {
                zpow <- c(zpow, rep(1, attr(zpow, "xnum")))
            }
            zpow <- cprtocmprpow(zpow[1], zpow[-1])
        }
        temp <- try(lookupvewa24ctqntl(
            prob = prob,
            ynum = 1, 
            zpow = zpow, 
            stat = "ct1",
            tolr = tolr
        ), TRUE)
        if (!inherits(temp, "try-error")) {
            if (!any(is.na(temp))) {
                return(drop(temp))
            }
        }
    }
	return(shinimcmprquant(prob, zpow, zadd, lgth, simu))
}