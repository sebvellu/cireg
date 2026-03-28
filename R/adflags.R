#' Select ADF lag length by information criteria
#'
#' Selects the number of lagged differences in the Augmented Dickey Fuller
#' regression using AIC or BIC. The function reuses the residuals produced
#' by adfstat to avoid rebuilding the regression repeatedly.
#' 
#' @details
#' The information criterion is computed as:
#'   smpl * log(rss / smpl) + penalty
#'
#' where rss is the sum of squared residuals from the ADF regression,
#' smpl is the effective sample size, and the penalty is:
#'   2 * k for AIC
#'   log(smpl) * k for BIC
#'
#' The parameter count k equals:
#'   1 for the lagged level term
#'   plus the number of lagged differences
#'   plus the number of deterministic regressors
#'
#' @param tsrs numeric vector. Time series to be tested for a unit root.
#' 
#' @param dvls optional numeric matrix. Deterministic regressors such as
#'   a constant or trend. Must have the same number of rows as `tsrs`.
#' 
#' @param maxl non negative integer. Maximum number of lagged differences
#'   considered. If `NULL`, a rule of thumb is used.
#' 
#' @param crit character string. Information criterion to use. Either
#' `"aic"` or `"bic"`.
#'
#' @return An integer. Selected number of lags. 
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tsrs <- cumsum(rnorm(200))
#' dvls <- cbind(rep(1, length(tsrs)))
#' adflags(tsrs, dvls, maxl = 12, crit = "bic")
#' }
#' 
#' @keywords internal

adflags <- function(tsrs, dvls = NULL, maxl = NULL, crit = "aic") {
	info <- Inf
	lvls <- NULL
	#
	if (is.null(maxl)) {
		maxl <- floor(12 * (length(tsrs)/100)^(1/4))
	}
	if (is.null(dvls)) {
		dnum <- 0
	} else {
        dvls <- as.matrix(dvls)
		dnum <- ncol(dvls)
	}
	#
    for (lags in 0:maxl) {
        rsds <- adfstat(tsrs, dvls, lags)$rsds
        smpl <- length(rsds)
        pars <- lags + 1 + dnum
        ssrr <- sum(rsds^2)
        if (crit == "aic") {
            temp <- 2 * pars
        } else { #if (crit == "bic") {
            temp <- log(smpl) * pars
        }
        temp <- smpl * log(ssrr/smpl) + temp
        if (temp < info) {
            lvls <- lags
            info <- temp
        }
    }
    #
    return(lvls)
}
