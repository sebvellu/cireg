getzvls <- function(zpow, xvls) {
    if (is.matrix(zpow)) {
        zvls <- cmprreg(zpow, xvls)
    } else {
        if (length(zpow) == 1) {
            zvls <- cbind(polydet(nrow(xvls), zpow), xvls)
        } else {
            # length(zpow) == ncol(xvls) + 1
            zpow <- cprtocmprpow(zpow[1], zpow[-1])
            zvls <- Recall(zpow, xvls)
        }
    }
    return(zvls)
}
