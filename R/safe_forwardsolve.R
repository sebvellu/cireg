#' Safe Forward Substitution
#'
#' Performs forward substitution on a lower triangular matrix.
#' The function first checks whether the diagonal elements of the
#' matrix are sufficiently large to allow stable use of `forwardsolve`.
#' If any diagonal element is smaller than a tolerance threshold,
#' forward substitution is considered unsafe and the function
#' automatically falls back to `qr_safe_solve`, which provides a
#' numerically stable solution using QR-based solving.
#'
#' This ensures that the function never fails due to singular or
#' near-singular triangular systems.
#'
#' @param lmat Lower triangular numeric matrix.
#' 
#' @param bmat Right-hand side vector or matrix.
#'
#' @return
#' A numeric vector or matrix containing the solution to the system.
#' The result is obtained either by `forwardsolve` or by `qr_safe_solve`,
#' depending on numerical safety.
#' 
#' @seealso qr_safe_solve
#'
#' @examples
#' lmat <- matrix(c(2, 0, 1, 3), 2, 2)
#' bmat <- c(1, 2)
#' safe_forwardsolve(lmat, bmat)
#' 
#' @keywords internal
#' 
safe_forwardsolve <- function(lmat, bmat) {
  tolr <- nrow(lmat) * .Machine$double.eps
  dvec <- diag(lmat)
  if (any(abs(dvec) < tolr)) {
    return(qr_safe_solve(lmat, bmat))
  }
 return(forwardsolve(lmat, bmat))
}
