#' Augmented Dickey Fuller test statistic (base implementation)
#'
#' Computes the Augmented Dickey Fuller t statistic for the coefficient
#' on the lagged level term using the Frisch Waugh Lovell theorem.
#' Lagged differences and deterministic regressors are partialled out
#' before estimating the coefficient on the lagged level.
#'
#' The function also returns the regression residuals from the full
#' ADF regression, which can be reused for information criteria based
#' lag selection.
#' 
#' @param tsrs numeric vector. Time series to be tested for a unit root.
#' 
#' @param dvls optional numeric matrix. Deterministic regressors such as
#' a constant or trend. Must have the same number of rows as `tsrs`.
#' 
#' @param lags non negative integer. Number of lagged differences to include.
#'
#' @return A list with components:
#' 
#'   - `stat`: numeric scalar. ADF t statistic for the lagged level term.
#'   - `rsds`: numeric vector. Regression residuals from the full ADF model.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tsrs <- cumsum(rnorm(200))
#' adfstat(tsrs, lags = 2)
#' }
#' 
#' @keywords internal

adfstat <- function(tsrs, dvls = NULL, lags = 0) {
	lgth <- length(tsrs)
	lmax <- lgth - 1
	ordr <- lags + 1
	#
	if (is.null(dvls)) {
		dvls <- matrix(0, nrow = lgth, ncol = 0)
	} else {
		dvls <- as.matrix(dvls)
	}
    #
	dtsr <- diff(tsrs)
	ltsr <- tsrs[-lgth]
	atsr <- matrix(0, lgth - ordr, lags)
	#
	for (indx in seq_len(lags)) {
		atsr[, indx] <- dtsr[(ordr - indx):(lmax - indx)]
	}
	#
	dtsr <- dtsr[ordr:lmax]
	ltsr <- ltsr[ordr:lmax]
	atsr <- cbind(atsr, dvls[(ordr + 1):lgth, , drop = FALSE])
	#
    if (ncol(atsr) > 0) {
	    qrat <- qr(atsr) # Apply Frish Waugh Theorem
	    ltsd <- qr.resid(qrat, ltsr)
	    dtsd <- qr.resid(qrat, dtsr)
    } else {
        ltsd <- ltsr
        dtsd <- dtsr
    }
	#
	deno <- sum(ltsd^2)
	coef <- sum(dtsd * ltsd)/deno
	rsds <- dtsd - coef * ltsd 
	shrt <- sum(rsds^2)/length(rsds) #lgth
	#
	return(list(stat = coef/sqrt(shrt/deno), rsds = rsds))
}
