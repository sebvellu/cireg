#' Fixed-b Wald-Type Tests for (Systems of) CMPRs Estimated by IM-OLS
#'
#' Performs a fixed-b Wald-type test based a result of a call to `imcmpr`.
#' 
#' @param sigl Significance level
#' 
#' @param objt Result of a call to `imcmpr`
#' 
#' @param rmat Restriction matrix
#' 
#' @param rvec Restriction vector
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

waldimcmprfbtest <- function(
    sigl, objt, rmat, rvec, smpl = 1000, simu = 10000,
    tolr = .Machine$double.eps
) {
    cffs <- c(t(objt$cffs))
    cvar <- objt$cvrb
    hypo <- nrow(rmat)
    rmat <- diag(1, objt$ynum) %x% rmat
    #
    stat <- waldstat(cffs, cvar, rmat, rvec)$stat
    crit <- getwfbimcmprqntl(
        prob = 1 - sigl,
        ynum = objt$ynum,
        zpow = objt$zpow,
        hypo = hypo,
        krnl = objt$krnl,
        bfrc = objt$fixb,
        zadd = objt$zadd,
        lgth = smpl,
        simu = simu,
        tolr = tolr
    )
    rjct <- c(stat > crit)
    return(list(stat = stat, crit = crit, rjct = rjct))
}
