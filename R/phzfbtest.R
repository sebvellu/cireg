#' Fixed-b Version of Phillips and Perron (1988) Z-Tests
#'
#' Performs the fixed-b version of the Phillips and Perrson (1988) Z-tests
#' provided by Vogelsang and Wagner (2013)
#' 
#' @param sigl Significance level.
#' 
#' @param tsrs Time series to be tested for a unit root.

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
#' @param step Integer in \{`1`, `2`\} specifying whether one- or two-step 
#' detrending should be applied. Either:
#' 
#'   - `1`: One-step detrending
#'   - `2`: Two-step detrending
#' 
#' Default is `2`. 
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
#' Vogelsang, T. J. and Wagner, M. (2013). A Fixed-b Perspective on the
#' Phillips-Perron Unit Root Tests. Econometric Theory 29, 609-628.
#' 
#' @export
#' 
phzfbtest <- function(
    sigl, tsrs, dpow = -1, krnl = "ba", band = "and", step = 2,
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    lgth <- length(tsrs)
    #bfrc <- 2 * 1:50/100
    dvls <- polydet(lgth, dpow)
    stat <- difustat(tsrs, dvls, NULL, step)
    stat <- phzfbstat(stat, krnl, band)
    #bfrc <- bfrc[which.min(abs(bfrc - stat[3]))]
    crit <- getphzfbqntl(sigl, dpow, step, krnl, stat[3], smpl, simu, tolr)
    stat <- stat[1:2]
    if (is.matrix(crit)) {
        rjct <- t(apply(crit, 1, function(x) {return(stat < x)}))
    } else {
        rjct <- (stat < crit)
    }
    return(list(stat = stat, crit = crit, rjct = rjct))
    #return(list(stat = stat, rjct = rjct))
}
