#' Phillips and Perron (1988) Z-Tests
#'
#' Performs Phillips and Perrson (1988) Z-tests
#' 
#' @param sigl Significance level.
#' 
#' @param tsrs Time series to be tested for a unit root.

#' @param dpow Integer specificing power of deterministics to consider (e.g.,
#' -1: no deterministics, 0: only constant; 1: constant and linear trend).
#' 
#' Default is -1.
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
#'   - `"and"`: Data-dependent rule (Andrews, 1991)
#'   - `"nw"`: Data-dependent rule (Newey & West, 1987)
#' 
#' Default is `"and"`.
#' 
#' @param step Integer, either 1 or 2. Specifys whether one- or two-step
#' detrending is applied.
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
#' Phillips, P. C. B. (1987). Time Series Regression with a Unit Root.
#' Econometrica 55, 277-301.
#' 
#' Phillips, P. C. B. and Perron, P. (1988). Testing for a Unit Root in Time
#' Series Regression. Biometrika 75, 335-346.
#' 
#' @export
#' 
phztest <- function(
    sigl, tsrs, dpow = -1, krnl = "ba", band = "and", step = 2,
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    lgth <- length(tsrs)
    dvls <- polydet(lgth, dpow)
    stat <- difustat(tsrs, dvls, NULL, step)
    stat <- phzstat(stat, krnl, band, step)
    crit <- getdifuqntl(sigl, dpow, 0, smpl, simu, tolr)
    if (is.matrix(crit)) {
        rjct <- t(apply(crit, 1, function(x) {return(stat < x)}))
    } else {
        rjct <- (stat < crit)
    }
    return(list(stat = stat, crit = crit, rjct = rjct))
    #return(list(stat = stat, rjct = rjct))
}
