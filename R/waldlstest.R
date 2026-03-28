#' Conventional Wald-Type Tests for (Systems of) CMPRs Estimated by
#' Integrated Modified Least Squares
#'
#' Performs a conventional Wald-type test based a result of a call to `imcmpr`.
#' 
#' @param sigl Significance level
#' 
#' @param objt Result of a call to `imls`
#' 
#' @param rmat Restriction matrix
#' 
#' @param rvec Restriction vector
#' 
#' @param rstr Logical. Should the Wald test be formed based for the free
#' parameters of the restricted estimator?
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

waldlstest <- function(sigl, objt, rmat, rvec, rstr = FALSE) {
    if (rstr) {
        cffs <- objt$cfff
        cvar <- objt$cvrf
    } else {
        cffs <- c(t(objt$cffs))
        cvar <- objt$cvrs
    }
    stat <- waldstat(cffs, cvar, rmat, rvec)
    crit <- stats::qchisq(1 - sigl, stat$degf)
    rjct <- c(stat$stat > crit)
    return(list(stat = stat$stat, crit = crit, rjct = rjct))
}
