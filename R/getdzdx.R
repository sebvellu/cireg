getdzdx <- function(zpow, xvls) {
    if (is.matrix(zpow)) {
	    dzdx <- cmprdzdx(xvls, zpow)
    } else {
        if (length(zpow) == 1) {
	        dzdx <- rbind(matrix(0, zpow + 1, ncol(xvls)), diag(1, ncol(xvls)))
	        dzdx <- replicate(nrow(xvls), dzdx)
        } else {
            zpow <- cprtocmprpow(zpow[1], zpow[-1])
	        dzdx <- Recall(zpow, xvls)
        }
    }
    return(dzdx)
}
