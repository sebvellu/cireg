 diffdtrnd <- function(xvls, dtrn = NULL, step = 2) {
	if (is.null(dtrn)) {
		vvls <- diff(xvls) #vdet
	} else {
		if (step == 1) {
			lgth <- nrow(xvls)
			qrdv <- qr(dtrn[-1, , drop = FALSE])
			vvls <- qr.resid(qrdv, xvls[-1, , drop = FALSE])
			vvls <- vvls - qr.resid(qrdv, xvls[-lgth, , drop = FALSE])
		} else { # if (step == 2) {
			vvls <- diff(qr.resid(qr(dtrn), xvls))
		}
	}
    return(vvls)
 }
