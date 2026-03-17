gettfbimcmprqntl <- function(
	prob, zpow, krnl, bfrc, zadd, 
    lgth = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    btbl <- 2 * 1:50/100
	btbl <- btbl[which.min(abs(btbl - bfrc))]
    if (is.null(zadd)) {
        if (!is.matrix(zpow)) {
            if (length(zpow) == 1) {
                temp <- try(lookupvowa14qntl(
                    prob = prob,
                    dpow = zpow, 
                    stat = "tsta",
                    xnum = attr(zpow, "xnum"),
                    hypo = 1,
                    krnl = krnl,
                    bfrc = btbl,
                    tolr = tolr
                ), TRUE)
                if (!inherits(temp, "try-error")) {
                    if (!any(is.na(temp))) {
                        return(temp)
                    }
                }
                zpow <- c(zpow, rep(1, attr(zpow, "xnum")))
            }
            zpow <- cprtocmprpow(zpow[1], zpow[-1])
        }
        temp <- try(lookupvewa24fbqntl(
            prob = prob,
            ynum = 1, 
            zpow = zpow, 
            stat = "tfb",
            hypo = 1,
            krnl = krnl,
            bfrc = btbl,
            tolr = tolr
        ), TRUE)
        if (!inherits(temp, "try-error")) {
            if (!any(is.na(temp))) {
                return(drop(temp))
            }
        }
    }
	return(tfbimcmprquant(
        prob = prob,
        zpow = zpow,
        krnl = krnl,
        bfrc = bfrc,
        zadd = zadd,
        lgth = lgth,
        simu = simu
    ))
}