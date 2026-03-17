#' Shin (1994) (D-OLS Based) Cointegration Test
#'
#' Performs a Shin (1994) Cointegration Test (CT), extended to the 
#' cointegrating polynomial regressions (CPR) setting (implicitly) by 
#' Knorre, Wagner and Grupe (2021).
#' 
#' @param sigl Significance level.
#' 
#' @param yvls Values of the integrated dependent variable.
#' 
#' @param zpow Specification of regression setting. Either:
#' 
#'   - An integer, in which case a linear 
#'     regression with highest polynomial time trend of specified power 
#'     is considered (e.g., -1: no deterministic components, 0: constant, 1: 
#'     constant and linear trend).
#'   - Vector of integers of size, in which case a 
#'     polynomial regression (PR) is considered with the first compoenent 
#'     specifying the hightest power of polynomial time trends
#'     and the remaining components specifying the highest powers of the
#'     non-deterministic regressors, i.e., columns of the matrix `xvls`.
#'   - A matrix, in which case a multivariate polynomial
#'     regression (MPR) is considered. Each row contains a multi-index 
#'     (i_0, i_1, ..., i_m) specifying the powers of the product 
#'     z_i = trend^(i_0) * x_1^(i_1) * ... * x_m^(i_m) considered as a 
#'     regressor in the regression. Here x_1, ..., x_m denote the columns
#'     of the matrix `xvls`. 
#'
#' @param xvls Matrix of integrated variables that do *not* cointegrate.
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
#' @param lale Information criterion to consider when choosing the 
#' optimal number of lags and leads for D-OLS estimation. 
#' Either `"aic"` or `"bic"`.
#' 
#' Default is `"aic"`.
#' 
#' @param dtrn Detrending of `xvls` to be applied in long-run covariance matrix 
#' estimation under conventional asymptotic theory. Either:
#' 
#'   - `NULL`: No detrening is applied
#'   - Matrix used as detrending regressors
#' 
#' Default is `NULL`.
#' 
#' @param step Integer in \{`1`, `2`\} specifying whether one- or two-step 
#' detrending in long-run covariance matrix estimation under conventional 
#' asymptotic theory is applied (only used, if `dtrn` is not `NULL`). Either:
#' 
#'   - `1`: One-step detrending
#'   - `2`: Two-step detrending
#' 
#' @param maxl Maximum number of lags and leads to consider for D-OLS
#' estimation. Either `NULL`, in which case an automated choice is made, an
#' integer specifying both the maximum number of lags and leads, 
#' or a vector of length two specifying the maximum number of lags 
#' and leads, respectively.
#' 
#' Default is `NULL`. 
#' 
#' @param symm Logical value indicating whether the number of leads
#' to choose should be equal to the number of lags.
#' 
#' Default is `TRUE`.
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
#' Knorre, F., Grupe, M. and Wagner, M. (2021). Monitoring Cointegrating
#' Polynomial Regressions: Theory and Application to the Environmental
#' Kuznets Curves for Carbon and Sulfur Dioxide Emissions. Econometrics 9, 12.
#' 
#' Shin, Y. (1994). A Residual-Based Test of the Null of Cointegration Against
#' the Alternative of No Cointegration. Econometric Theory 10, 91-115.
#'  
#' @export
  
shincmprtest <- function(
    sigl, yvls, zpow, xvls, krnl = "ba", band = "and", lale = "aic",
    dtrn = NULL, step = 2, maxl = NULL, symm = TRUE,
    smpl = 1000, simu = 10000, tolr = .Machine$double.eps
) {
    lgth <- length(yvls)
    attr(zpow, "xnum") <- ncol(xvls)
    stat <- shincmprstat(
        yvls = yvls,
        zpow = zpow,
        xvls = xvls,
        lale = lale,
        krnl = krnl,
        band = band,
        dtrn = dtrn,
        step = step,
        maxl = maxl,
        symm = symm
    )
    crit <- getshincmprqntl(1 - sigl, zpow, smpl, simu, tolr)
    rjct <- (stat > crit)
    return(list(stat = stat, crit = crit, rjct = rjct))
}
