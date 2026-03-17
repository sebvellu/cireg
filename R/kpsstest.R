#' Kwiatkowski, Phillips, Schmidt and Shin (1992) Test
#'
#' Performs a KPSS test.
#' 
#' @param sigl Significance level.
#' 
#' @param tsrs Time series to be tested for stationarity.

#' @param dpow Integer specificing power of deterministics to consider (e.g.,
#' `-1`: no deterministics, `0`: only constant; `1`: constant and linear trend).
#' 
#' Default is `-1`.
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
#' 
#' Default is `"and"`.
#' 
#' @param smpl Sample size to consider in case critical values 
#' have to be simulated.
#' 
#' @param simu Number of repetitions in case critical values have to 
#' be simulated.
#' 
#' @param tolr Tolerance limit for lookup tables of critical values.
#'
#' @return A list containing:
#' 
#'   - `stat`: Test statistic
#'   - `crit`: Critical value
#'   - `rjct`: Logical value indicating the rejection of null hypothesis
#' 
#' @references
#' Kwiatkowski, D., Phillips, P. C. B., Schmidt, P. and Shin, Y. (1992).
#' Testing the Null Hypothesis of Stationarity Against the Alternative of a
#' Unit Root. Journal of Econometrics 54, 159-178.
#' 
#' @export
#' 
kpsstest <- function(
    sigl, tsrs, dpow = NULL, krnl = "ba", band = "and",
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    lgth <- length(tsrs)
    dvls <- polydet(lgth, dpow)
    stat <- kpssstat(tsrs, dvls, krnl, band)
    crit <- getkpssqntl(1 - sigl, dpow, smpl, simu, tolr)
    rjct <- (stat > crit)
    return(list(stat = stat, crit = crit, rjct = rjct))
}
