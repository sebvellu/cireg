#' Safe Backward Substitution
#'
#' Performs backward substitution on an upper triangular matrix.
#' The function first checks whether the diagonal elements of the
#' matrix are sufficiently large to allow stable use of `backsolve`.
#' If any diagonal element is smaller than a tolerance threshold,
#' backward substitution is considered unsafe and the function
#' automatically falls back to `qr_safe_solve`, which provides a
#' numerically stable solution using QR-based solving.
#'
#' This ensures that the function never fails due to singular or
#' near-singular triangular systems.
#'
#' @param rmat Upper triangular numeric matrix.
#' 
#' @param bmat Right-hand side vector or matrix.
#'
#' @return
#' A numeric vector or matrix containing the solution to the system.
#' The result is obtained either by `backsolve` or by `qr_safe_solve`,
#' depending on numerical safety.
#' 
#' @seealso qr_safe_solve
#'
#' @examples
#' \dontrun{
#' rmat <- matrix(c(2, 1, 0, 3), 2, 2)
#' bmat <- c(1, 2)
#' safe_backsolve(rmat, bmat)
#' }
#' 
#' @keywords internal
#' 
safe_backsolve <- function(rmat, bmat) {
  tolr <- nrow(rmat) * .Machine$double.eps
  dvec <- diag(rmat)
  if (any(abs(dvec) < tolr)) {
    return(qr_safe_solve(rmat, bmat))
  }
  return(backsolve(rmat, bmat))
}
