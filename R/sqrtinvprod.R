sqrtinvprod <- function(qrwv, foba = TRUE) { # let qrwv be qr decomp of wvls then estimates 
    if (foba) {
        rqrw <- qr.R(qrwv)
        #whiv <- qr_safe_solve(rqrw, qr_safe_solve(t(rqrw), diag(1, ncol(rqrw))))
        temp <- safe_forwardsolve(t(rqrw), diag(1, ncol(rqrw)))
        whiv <- safe_backsolve(rqrw, temp)
    } else {
	    rqrw <- qr.R(qrwv)
	    #whiv <- qr_safe_solve(t(rqrw), diag(1, ncol(rqrw)))
	    whiv <- safe_forwardsolve(t(rqrw), diag(1, ncol(rqrw)))
    }
    return(whiv)
}