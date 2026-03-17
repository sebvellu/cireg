tstat <- function(cffs, cvar, rvec, rval) {
    tsta <- sum(rvec * c(cffs)) - rval
    midl <- sum(rvec * (cvar %*% rvec))
    tsta <- tsta/sqrt(midl)
    return(tsta)
}