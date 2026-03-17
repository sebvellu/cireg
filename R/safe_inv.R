#' Safely invert a matrix
#'
#' Tries to invert a matrix using the solve function. If `solve` fails,
#' usually because the matrix is singular or ill-conditioned,
#' the function returns the result of `mp_inv` instead.
#'
#' @param matr A numeric matrix to invert.
#'
#' @return The inverse of matr, or its Moore-Penrose pseudoinverse
#' if solve fails.
#'
#' @examples
#' safe_inv(diag(3))
#'
#' @seealso mp_inv
#' 
#' @keywords internal
#' 
safe_inv <- function(matr) {
	tsol <- try(solve(matr), TRUE)
	if (inherits(tsol, "try-error")) {
		return(mp_inv(matr))
	} else {
		return(tsol)
	}
}