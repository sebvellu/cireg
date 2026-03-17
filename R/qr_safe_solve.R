#' Safe Linear Solver Using QR or Pseudoinverse Fallback
#'
#' This function solves the linear system A x = B using `qr.solve`.
#' If `qr.solve` fails (for example due to rank deficiency), the function
#' falls back to using the Moore-Penrose pseudoinverse computed by
#' `mp_inv`.
#'
#' @details
#' The function first attempts:
#' 
#' `qr.solve(mata, matb)`
#' 
#' If this produces an error, a pseudoinverse-based solution is used:
#' 
#' `mp_inv(mata) %*% matb`
#'
#' @param mata Numeric matrix A in the linear system A x = B.
#' 
#' @param matb Numeric matrix or vector B in the linear system A x = B.
#'
#' @return
#' A numeric matrix or vector giving the solution to A x = B. If
#' `qr.solve` succeeds, its result is returned. If it fails, the solution
#' computed using `mp_inv(A) %*% B` is returned.
#' 
#' @seealso mp_inv
#'
#' @examples
#' A <- matrix(c(1, 2, 2, 4), 2)
#' b <- c(1, 1)
#' qr_safe_solve(A, b)
#' 
#' @keywords internal
#' 
qr_safe_solve <- function(mata, matb) {
    tsol <- try(qr.solve(mata, matb), TRUE)
    if (inherits(tsol, "try-error")) {
        return(mp_inv(mata) %*% matb)
    } else {
        return(tsol)
    }
}
