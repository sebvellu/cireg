getctimcmprqntl <- function(
	prob, ynum, zpow, zadd,
    lgth = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    if (is.null(zadd)) {
        tmpa <- try(lookupvewa24ctqntl(
            prob = prob,
            ynum = ynum, 
            zpow = zpow, 
            stat = "ct1",
            tolr = tolr
        ), TRUE)
        if (!inherits(tmpa, "try-error")) {
            if (!any(is.na(tmpa))) {
                tmpb <- try(lookupvewa24ctqntl(
                    prob = prob,
                    ynum = ynum, 
                    zpow = zpow, 
                    stat = "ct2",
                    tolr = tolr
                ), TRUE)
                tmpc <- try(lookupvewa24ctqntl(
                    prob = prob,
                    ynum = ynum, 
                    zpow = zpow, 
                    stat = "ct3",
                    tolr = tolr
                ), TRUE)
                #
                return(drop(cbind(tmpa, tmpb, tmpc)))
            }
        }
    }
	return(ctimcmprquant(prob, ynum, zpow, zadd, lgth, simu))
}