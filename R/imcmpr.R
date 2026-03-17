#' Integrated Modified Estimation of Systems of Cointegrating Multivariate
#' Polynomial Regressions (IM-SCMPR)
#'
#' Computes IM-OLS and IM-GLS estimates for a system of
#' cointegrating multivariate polynomial regressions.
#'
#' @details
#' Model setup:
#' 
#' y = cf_mat * Z + u
#' 
#' Z = \[z_1, ..., z_I\]'
#' 
#' y = \[y_1, ..., y_n\]'
#' 
#' with: 
#' 
#' z_i = trend^(i_0) * x_1^(i_1) * ... * x_m^(i_m) \cr
#' 
#' IM (augmented partial sum) regression considered:
#' 
#' S^y = cf_mat * S^Z + cf_aug * x + S^u
#' 
#' x = \[x_1, ..., x_m\]'
#' 
#' where:
#' 
#' vec(\[cf_mat, cf_aug\]') = D * cf_free + d \cr
#' 
#' If `zadd` is not `NULL`, the IM regression considered becomes:
#' 
#' S^y = cf_mat * S^Z + cf_aug * \[w', x'\]' + S^u
#' 
#' with w corresponding to the regressors specified by `zadd`.
#' 
#' @param yvls Matrix of size T x n of integrated variables that cointegrate 
#' with `xvls` (see y in the *Details* section); left-hand side of the system 
#' (dependent variables). 
#' 
#' @param zpow Specification of regression setting. Either:
#' 
#'   - An integer, in which case a linear cointegrating 
#'     regression with highest polynomial time trend of specified power 
#'     is considered (e.g., -1: no deterministic components, 0: constant, 1: 
#'     constant and linear trend).
#'   - Vector of integers of size (m + 1), in which case a cointegrating 
#'     polynomial regression (CPR) is considered with the first compoenent 
#'     specifying the hightest power of polynomial time trends
#'     and the remaining components specifying the highest powers of the
#'     non-deterministic regressors, i.e., columns of the matrix `xvls`.
#'   - Matrix of size I x (m + 1), in which case a cointegrating multivariate
#'     polynomial regression (CMPR) is considered. The i-th row contains
#'     the multi-index (i_0, i_1, ..., i_m) specifying the powers of the product 
#'     z_i = trend^(i_0) * x_1^(i_1) * ... * x_m^(i_m) considered as a 
#'     regressor in the regression.
#' 
#' @param xvls Matrix of size T x m of integrated variables that do *not*
#' cointegrate (see x in the *Details* section).
#' 
#' @param zadd Matrix of size T x q of regressors to be added to the IM 
#' (augmented partial sum) regression (see w in the *Details* section).
#' Either:
#'   
#'   - `NULL`: Empty matrix of size T x 0 (no additional regressors)
#'   - Matrix of size T x q
#' 
#' Default is `NULL`.
#' 
#' @param rmat Restriction matrix of size ((I + m) * n) x g, denoted as D in 
#' the *Details* section, such that vec(\[cf_mat, cf_aug\]') = D * cf_free + d.
#' Either:
#' 
#'   - `NULL`: Identity matrix of size ((I + m) * n) x ((I + m) * n)
#'   - Matrix of size ((I + m) * n) x g 
#' 
#' Default is `NULL`.
#' 
#' @param rvec Restriction vector of size (I + m) * n, denoted as d in the 
#' *Details* section, such that vec(\[cf_mat, cf_aug\]') = D * cf_free + d.
#' Either:
#' 
#'   - `NULL`: Zero vector of size (I + m) * n
#'   - Vector of size (I + m) * n
#' 
#' Default is `NULL`.
#' 
#' @param wght Weighting matrix for IM-GLS estimation. Either:
#' 
#'   - `NULL`: Identity matrix of size n x n
#'   - `"PaOg1"`: Use estimate of Omega_uu^(-1)
#'   - `"PaOg2"`: Use estimate of Omega_u.v^(-1)
#'   - Matrix of size n x n
#' 
#' Default is `NULL`.
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
#'   - `"nwt"`: Rule-of-thumb, i.e., floor(4 * (T / 100)^(2 / 9))
#' 
#' Default is `"and"`.
#' 
#' @param dtrn Detrending of `xvls` to be applied in long-run covariance matrix 
#' estimation under conventional asymptotic theory. Either:
#' 
#'   - `NULL`: No detrening is applied
#'   - Matrix of size T x p used as detrending regressors
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
#' Default is `2`. 
#'
#' @return A list containing:
#' 
#'   - `cffs`: Matrix of size n x (I + m) of unrestricted IM-OLS estimates of 
#'     \[cf_mat, cf_aug\].
#'   - `fits`: Matrix of size T x n of fitted values from unrestricted IM-OLS 
#'     estimation.
#'   - `rsds`: Matrix of size T x n of residuals from unrestricted IM-OLS 
#'     estimation.
#'   - `clrv`: Matrix of size n x n of the estimated conditional long-run
#'     covariance matrix Omega_u.v used for conventional asymptotic theory.
#'   - `cvrs`: Matrix of size ((I + m) * n) x ((I + m) * n) giving the 
#'     estimated asymptotic covariance matrix of vec(`cffs`') for conventional
#'     asymptotic theory.
#'   - `csds`: Matrix of size n x (I + m) giving the estimated asymptotic
#'     standard deviations of `cffs` for conventional asymptotic theory.
#'   - `rsdm`: Matrix of size T x n of modified IM-OLS residuals used for
#'     fixed-b long-run covariance matrix estimation.
#'   - `cvrb`: Matrix of size ((I + m) * n) x ((I + m) * n) giving the 
#'     estimated asymptotic covariance matrix of vec(`cffs`') for fixed-b
#'     asymptotic theory.
#'   - `csdb`: Matrix of size n x (I + m) giving the estimated asymptotic 
#'     standard deviations of `cffs` for fixed-b asymptotic theory.
#'   - `fixb`: Bandwidth-to-sample-size ratio. If `band` is  either `and` or
#'     `nw`, computed usind the first differences of `rsdm`. Otherwise,
#'     if `band` is either `nwt` or an integer in \{1, ..., T\},
#'     the corresponding bandwidth divided by the sample size.
#'   - `cffr`: Matrix of size n x (I + m) giving the restricted IM-GLS
#'     estimates of \[cf_mat, cf_aug\] with 
#'     vec(\[cf_mat, cf_aug\]') = D * cf_free + d.
#'   - `fitr`: Matrix of size T x n of fitted values from restricted IM-GLS 
#'     estimation.
#'   - `rsdr`: Matrix of size T x n of residuals from restricted IM-GLS 
#'     estimation.
#'   - `cvrr`: Matrix of size ((I + m) * n) x ((I + m) * n) giving the 
#'     estimated asymptotic covariance matrix of vec(`cffr`') for conventional 
#'     asymptotic theory.
#'   - `csdr`: Matrix of size n x (I + m) giving the estimated asymptotic
#'     standard deviations of `cffr` for conventional asymptotic theory.
#'   - `cfff`: Vector of size g giving the IM coefficient estimates of cf_free.
#'   - `cvrf`: Matrix of size g x g giving the estimated asymptotic covariance 
#'     matrix of `cfff` for conventional asymptotic theory.
#'
#' @note
#' The function does not verify asymptotic rank conditions or check that
#' the restrictions D and d are correctly specified.
#' It assumes that the restrictions do not involve elements of cf_aug.
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

imcmpr <- function(
	yvls, zpow, xvls, zadd = NULL, rmat = NULL, rvec = NULL, wght = NULL, 
	krnl = "ba", band = "and", dtrn = NULL, step = 2
) {
	zvls <- getzvls(zpow, xvls)
	rslt <- imls(
		yvls = yvls,
		zvls = zvls,
		xvls = xvls,
		zadd = zadd,
		rmat = rmat,
		rvec = rvec,
		wght = wght,
		krnl = krnl,
		band = band,
		dtrn = dtrn,
		step = step
	)
	rslt$zpow <- zpow
    attr(rslt$zpow, "xnum") <- ncol(xvls)
	rslt$ynum <- NCOL(yvls)
	rslt$krnl <- krnl
	rslt$zadd <- zadd
	return(rslt)
}