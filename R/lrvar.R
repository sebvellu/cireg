#' Estimation of the Long-Run, Half Long-Run and One-Sided Long-Run 
#' Covariance Matrix of a Linear process
#' 
#' Computes long-run, half long-run and one-sided long-run covariance matrix
#' estimates of a (multivariate) time series
#'
#' @param tsrs A vector or matrix of values for the time series.
#' 
#' @param krnl Kernel function used for long-run variance estimation. Either:
#' 
#'   - `"tr"`: Truncated
#'   - `"ba"`: Bartlett
#'   - `"pa"`: Parzen
#'   - `"bo"`: Bohman
#'   - `"da"`: Daniell
#'   - `"qs"`: Quadratic Spectral
#' 
#' Default is `"ba"`.
#' 
#' @param band Bandwidth specification. Either:
#' 
#'   - Integer in \{1, ..., T\} (explicit bandwidth value)
#'   - `"and"`: Data-dependent rule (Andrews, 1991)
#'   - `"nw"`: Data-dependent rule (Newey & West, 1987)
#'   - `"nwt"`: Rule-of-thumb, i.e., floor(4 * (T / 100)^(2 / 9))
#' 
#' Default is `"and"`.
#' 
#' @return A list containing:
#' 
#'   - `bwdh`: The bandwidth computed/specified
#'   - `halfvar`: The half long-run covariance matrix
#'   - `shrtvar`: The contemporaenous covariance matrix
#'   - `longvar`: The long-run covariance matrix
#'   - `deltvar`: The one-sided long-run covariance matrix, i.e.,
#'   `halfvar` + `shrtvar`
#'
#' @seealso bwclc, lrwghts, computelrvec, computelrmat
#' 
#' @references
#' Andrews, D. W. K. (1991).
#' "Heteroskedasticity and autocorrelation consistent covariance matrix 
#' estimation." *Econometrica*, 59(3), 817-858.
#'
#' Newey, W. K., and West, K. D. (1987).
#' "A simple, positive semi-definite, heteroskedasticity and autocorrelation 
#' consistent covariance matrix." *Econometrica*, 55(3), 703-708.
#' 
#' @keywords internal
#' 
lrvar <- function(tsrs, krnl = "ba", band = "and") {
	lgth <- NROW(tsrs)
	nums <- NCOL(tsrs)
	#
	bwdh <- bwclc(tsrs, krnl, band, NULL)
	wght <- lrwghts(lgth, bwdh, krnl)
	#
	if (nums == 1) {
		lvar <- computelrvec(tsrs/sqrt(lgth), wght) #computelrvecC
	} else {
		lvar <- computelrmat(tsrs/sqrt(lgth), wght)
	}
	#
	return(list(
		bwdh = bwdh,
		shrtvar = lvar$shrt,
		halfvar = lvar$half,
		longvar = lvar$long,
		deltvar = lvar$delt
	))
}
