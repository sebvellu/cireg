#' CT Cointegration Tests for (Systems of) CMPRs Estimated by IM-OLS
#'
#' Performs CT tests for (systems of) CMPRs based a result of a call 
#' to `imcmpr`.
#' 
#' @param sigl Significance level
#' 
#' @param objt Result of a call to `imcmpr`
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
#' Nyblom, J. and Harvey, A. (2000). Tests of Common Stochastic Trends.
#' Econometric Theory 16, 176-199.
#' 
#' Shin, Y. (1994). A Residual-Based Test of the Null of Cointegration Against
#' the Alternative of No Cointegration. Econometric Theory 10, 91-115.
#' 
#' Veldhuis, S. and Wagner, M. (2026). Estimation and Inference for
#' Cointegrated Systems of Multi-Factor Production Functions: Modelling the
#' Joint Behavior of GDP and Emissions. Mimeo.
#' 
#' @export
#' 
ctimcmprtest <- function(
	sigl, objt, smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    crsd <- crossprod(objt$rsds)
    ccrs <- crossprod(apply(objt$rsds, 2, cumsum))
    stat <- ctstatint(objt$rsds, objt$clrv)
    crit <- getctimcmprqntl(1 - sigl, objt$ynum, objt$zpow, objt$zadd, smpl, simu, tolr)
    if (is.matrix(crit)) {
        rjct <- t(apply(crit, 1, function(x) {return(stat > x)}))
    } else {
        rjct <- (stat > crit)
    }
    return(list(stat = stat, crit = crit, rjct = rjct))
}
