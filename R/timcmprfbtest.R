#' Fixed-b t-Type Tests for CMPRs Estimated by IM-OLS
#'
#' Performs a fixed-b t-type test based a result of a call to `imcmpr`.
#' 
#' @param sigl Significance level
#' 
#' @param objt Result of a call to `imcmpr` with `ncol(yvls) == 1`
#' 
#' @param rvec Restriction vector
#' 
#' @param rval Restriction value
#' 
#' @param dirc Integer specifying the test direction: -1 for left-tailed,
#' 1 for right-tailed, or 0 for two-sided.
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

timcmprfbtest <- function(
    sigl, objt, rvec, rval, dirc = 0, smpl = 1000, simu = 10000,
    tolr = .Machine$double.eps
) {
    stopifnot(objt$ynum == 1)
    #
    cffs <- c(t(objt$cffs))
    cvar <- objt$cvrb
    #
    stat <- tstat(cffs, cvar, rvec, rval)
    cfun <- function(prob) {
        rslt <- gettfbimcmprqntl(
            prob = prob,
            zpow = objt$zpow,
            krnl = objt$krnl,
            bfrc = objt$fixb,
            zadd = objt$zadd,
            lgth = smpl,
            simu = simu,
            tolr = tolr
        )
        return(rslt)
    }
    if (dirc == 0) {
        crit <- cfun(1 - sigl/2)
        rjct <- c(abs(stat) > crit)
    } else if (dirc == 1) {
        crit <- cfun(1 - sigl)
        rjct <- c(stat > crit)
    } else {
        crit <- -cfun(1 - sigl)
        rjct <- c(stat < crit)
    }
    return(list(stat = stat, crit = crit, rjct = rjct))
}
