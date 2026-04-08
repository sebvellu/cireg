#' Augmented Dickey and Fuller (1979) Test
#'
#' Performs an augmented Dickey Fuller (1979) (ADF) Test
#' 
#' @param sigl Significance level.
#' 
#' @param tsrs Time series to be tested for a unit root.

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
#' @export
#' 
adftest <- function(
    sigl, tsrs, dpow = -1, maxl = NULL, crit = "aic",
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    lgth <- length(tsrs)
    dvls <- helperkit::polydet(lgth, dpow)
    lags <- adflags(tsrs, dvls, maxl, crit)
    stat <- adfstat(tsrs, dvls, lags)$stat
    crit <- getdifuqntl(sigl, dpow, 0, smpl, simu, tolr)
    if (is.matrix(crit)) {
        crit <- crit[, 2]
    } else {
        crit <- crit[2]
    }
    rjct <- (stat < crit)
    return(list(stat = stat, crit = crit, rjct = rjct))
    #return(list(stat = stat, rjct = rjct))
}
