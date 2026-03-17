calcbhlf <- function(cvls, clrv, rmat, wght) {
    bhlf <- ((safe_chol(clrv) %*% wght) %x% cvls) %*% rmat
    return(bhlf)
}
