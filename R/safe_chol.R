#' Safe Cholesky Decomposition
#'
#' Computes a Cholesky-like factor of a variance-covariance matrix.
#' The function first attempts to compute the standard Cholesky
#' decomposition using `chol`. If this succeeds, the result is returned.
#' If `chol` fails (for example because the matrix is ill-conditioned
#' or not numerically positive definite), the function falls back to an
#' eigenvalue-based method. Eigenvalues are clipped to a minimum value
#' to ensure numerical stability, and a substitute factor is constructed.
#'
#' This function always returns a real matrix A such that t(A) %*% A
#' is a positive semidefinite approximation of the input matrix.
#'
#' @param matr A numeric square matrix, typically a covariance matrix.
#'
#' @return
#' A numeric matrix. If `chol` succeeds, this is the upper triangular
#' Cholesky factor. Otherwise it is a stable eigenvalue-based factor.
#'
#' @examples
#' # Basic usage
#' matr <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
#' safe_chol(matr)
#'
#' # Ill-conditioned matrix example
#' matr <- matr
#' matr[1, 1] <- 1e-14
#' safe_chol(matr)
#' 
#' @keywords internal
#' 
safe_chol <- function(matr) {
    tolr <- nrow(matr) * .Machine$double.eps
    matr <- (matr + t(matr))/2 # enforce symmetry
    # Try standard chol
    rslt <- try(chol(matr), silent = TRUE)
    if (!inherits(rslt, "try-error")) {
  	    return(rslt)
    }
    # Fallback: eigenvalue clipping
    eign <- eigen(matr, symmetric = TRUE)
    vals <- pmax(eign$values, tolr)
    rslt <- diag(sqrt(vals), nrow(matr)) %*% t(eign$vectors)
    return(rslt)
}
