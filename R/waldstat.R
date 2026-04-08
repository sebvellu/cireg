waldstat <- function(cffs, cvar, rmat, rvec) {
    wald <- rmat %*% cffs - rvec
    midl <- helperkit::safesolve(tcrossprod(rmat %*% cvar, rmat))
    wald <- drop(crossprod(wald, midl %*% wald))
    wadf <- length(rvec)
    return(list(stat = wald, degf = wadf))
}
