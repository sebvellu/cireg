getwfbimcmprqntl <- function(
	prob, ynum, zpow, hypo, krnl, bfrc, zadd, 
    lgth = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    btbl <- 2 * 1:50/100
	btbl <- btbl[which.min(abs(btbl - bfrc))]
    if (is.null(zadd)) {
        if (!is.matrix(zpow)) {
            if (length(zpow) == 1) {
                if (ynum == 1) {
                    temp <- try(lookupvowa14qntl(
                        prob = prob,
                        dpow = zpow,
                        stat = "wald",
                        xnum = attr(zpow, "xnum"),
                        hypo = hypo,
                        krnl = krnl,
                        bfrc = btbl,
                        tolr = tolr
                    ), TRUE)
                    if (!inherits(temp, "try-error")) {
                        if (!any(is.na(temp))) {
                            return(temp)
                        }
                    }
                }
                zpow <- c(zpow, rep(1, attr(zpow, "xnum")))
            }
            zpow <- cprtocmprpow(zpow[1], zpow[-1])
        }
        if (ynum == 1) {
            temp <- try(lookupvowa24qntl(
                prob = prob,
                zpow = zpow, 
                hypo = hypo,
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
        temp <- try(lookupvewa24fbqntl(
            prob = prob,
            ynum = ynum, 
            zpow = zpow, 
            stat = "wfb",
            hypo = hypo,
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
	return(wfbimcmprquant(
        prob = prob, 
        ynum = ynum,
        zpow = zpow,
        hypo = hypo,
        krnl = krnl,
        bfrc = bfrc,
        zadd = zadd,
        lgth = lgth,
        simu = simu
    ))
}