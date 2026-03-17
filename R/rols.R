 rols <- function(yvls, xvls, spec = "ols") { # restricted ols
	yvls <- as.matrix(yvls)
	xvls <- as.matrix(xvls)
	#
	lgth <- nrow(yvls)
	ynum <- ncol(yvls)
	xnum <- ncol(xvls)
	#
	qrxv <- qr(xvls)
	xvri <- backsolid(qr.R(qrxv))
	cffs <- crossprod(yvls, tcrossprod(qr.Q(qrxv), xvri)) # == crossprod(yvls, xvls) %*% solve(crossprod(xvls))
	rsds <- yvls - tcrossprod(xvls, cffs)
	#
	if (identical(spec, "ols")) {
		return(list(cffs = cffs, rsds = rsds))
	}
	#
	rmat <- spec[[1]] # restrictions on vec(Phi')
	rvec <- spec[[2]]
	wght <- spec[[3]]
	whlf <- chol(wght)
	#
	temp <- qr.R(qrxv) %x% whlf
	qrtp <- qr(temp %*% rmat)
	tpri <- backsolid(qr.R(qrtp))
	temq <- as.vector(whlf %*% crossprod(yvls, qr.Q(qrxv)))
	cfre <- tpri %*% crossprod(qr.Q(qrtp), temq - temp %*% rvec)
	#
	cffr <- rmat %*% cfre + rvec
	cffr <- matrix(cffr, ynum, xnum)
	rsdr <- yvls - tcrossprod(xvls, cffr)
	#
	return(list(cfre = cfre, cffr = cffr, rsdr = rsdr))
}