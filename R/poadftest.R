#' Phillips-Ouliaris (1990) Version of ADF Test
#'
#' Performs a non-cointegration Phillips-Ouliaris (1990) ADF test.
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
#' @param maxl Maximum number of lags to consider. Either `NULL`
#' in which case an automated choice is made, or an integer 
#' specifying the maximum number of lags.
#' 
#' Default is `NULL`. 
#' 
#' @param crit Information criterion to consider when choosing the 
#' optimal number of lags. Either `"aic"` or `"bic"`.
#' 
#' Default is `"aic"`.
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
#' Dickey, D. A. and Fuller, W. A. (1979). Distribution of the Estimators for
#' Autoregressive Time Series With a Unit Root. Journal of the American
#' Statistical Association 74, 427-431.
#' 
#' Phillips, P. C. B. and Ouliaris, S. (1990). Asymptotic Properties of
#' Residual Based Tests for Cointegration. Econometrica 58, 165-193.
#' 
#' @export

poadftest <- function(
    sigl, yvls, xvls, dpow = -1, maxl = NULL, crit = "aic",
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    lgth <- length(yvls)
    dvls <- polydet(lgth, dpow)
    tsrs <- qr.resid(qr(cbind(dvls, xvls)), yvls)
    lags <- adflags(tsrs, NULL, maxl, crit)
    stat <- adfstat(tsrs, NULL, lags)$stat
    crit <- getdifuqntl(sigl, dpow, NCOL(xvls), smpl, simu, tolr)
    if (is.matrix(crit)) {
        crit <- crit[, 2]
    } else {
        crit <- crit[2]
    }
    rjct <- (stat < crit)
    return(list(stat = stat, crit = crit, rjct = rjct))
    #return(list(stat = stat, rjct = rjct))
}
