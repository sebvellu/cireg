#' Phillips-Ouliaris (1990) Version of Phillips and Perron (1988) Z-Tests
#'
#' Performs non-cointegration Phillips-Ouliaris (1990) Z-tests
#' 
#' @param sigl Significance level.
#' 
#' @param yvls Values of the integrated dependent variable.
#' 
#' @param xvls Matrix of integrated regressors that do *not* cointegrate.
#' 
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
#'   - `stat`: Test statistics
#'   - `crit`: Critical values
#'   - `rjct`: Logical values indicating the rejection of null hypothesis
#' 
#' @references 
#' Phillips, P. C. B. and Ouliaris, S. (1990). Asymptotic Properties of
#' Residual Based Tests for Cointegration. Econometrica 58, 165-193.
#' 
#' Phillips, P. C. B. and Perron, P. (1988). Testing for a Unit Root in Time
#' Series Regression. Biometrika 75, 335-346.
#' 
#' @export
#' 
poztest <- function(
    sigl, yvls, xvls, dpow = -1, krnl = "ba", band = "and",
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    lgth <- length(yvls)
    dvls <- polydet(lgth, dpow)
    stat <- difustat(yvls, dvls, xvls, NULL)
    stat <- phzstat(stat, krnl, band, NULL)
    crit <- getdifuqntl(sigl, dpow, NCOL(xvls), smpl, simu, tolr)
    if (is.matrix(crit)) {
        rjct <- t(apply(crit, 1, function(x) {return(stat < x)}))
    } else {
        rjct <- (stat < crit)
    }
    return(list(stat = stat, crit = crit, rjct = rjct))
    #return(list(stat = stat, rjct = rjct))
}
