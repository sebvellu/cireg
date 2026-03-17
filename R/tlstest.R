#' Conventional t-Type Tests for CMPRs Estimated by
#' Integrated Modified Least Squares
#'
#' Performs a conventional t-type test based a result of a call to `imcmpr`.
#' 
#' @param sigl Significance level
#' 
#' @param objt Result of a call to `imcmpr`
#' 
#' @param rvec Restriction vector
#' 
#' @param rval Restriction value
#' 
#' @param rstr Logical. Should the Wald test be formed based for the free
#' parameters of the restricted estimator?
#' 
#' @param dirc Integer specifying the test direction: -1 for left-tailed,
#' 1 for right-tailed, or 0 for two-sided.
#'
#' @return A list containing:
#' 
#'   - `stat`: Test statistic
#'   - `crit`: Critical value
#'   - `rjct`: Logical value indicating the rejection of null hypothesis
#' 
#' @references
#' Veldhuis, S. and Wagner, M. (2026). Estimation and Inference for
#' Cointegrated Systems of Multi-Factor Production Functions: Modelling the
#' Joint Behavior of GDP and Emissions. Mimeo.
#' 
#' Vogelsang, T. J. and Wagner, M. (2014). Integrated Modified OLS Estimation
#' and Fixed-b Inference for Cointegrated Regressions. Journal of
#' Econometrics 178, 741-760.
#' 
#' Vogelsang, T. J. and Wagner, M. (2024). Integrated Modified OLS Estimation
#' and Fixed-b Inference for Cointegrating Multivariate Polynomial
#' Regressions. IHS Working Paper Series 53.
#' 
#' @export

tlstest <- function(sigl, objt, rvec, rval, rstr = FALSE, dirc = 0) {
    if (rstr) {
        cffs <- objt$cfff
        cvar <- objt$cvrf
    } else {
        cffs <- c(t(objt$cffs))
        cvar <- objt$cvrs
    }
    stat <- tstat(cffs, cvar, rvec, rval)
    if (dirc == 0) {
        crit <- qnorm(1 - sigl/2)
        rjct <- c(abs(stat) > crit)
    } else if (dirc == 1) {
        crit <- qnorm(1 - sigl)
        rjct <- c(stat > crit)
    } else {
        crit <- -qnorm(1 - sigl)
        rjct <- c(stat < crit)
    }
    return(list(stat = stat, crit = crit, rjct = rjct))
}
