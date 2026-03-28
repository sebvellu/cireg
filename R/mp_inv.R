#' Moore-Penrose Pseudoinverse
#'
#' This function computes the Moore-Penrose pseudoinverse of a numeric
#' matrix using its singular value decomposition (SVD). Very small
#' singular values are treated as zero using a tolerance based on
#' machine precision.
#'
#' @details
#' The matrix is decomposed as U D V' using svd. Singular values smaller
#' than
#' 
#' `max(dim(matr)) * .Machine$double.eps * max(svd$d)`
#' 
#' are discarded. The pseudoinverse is constructed using the remaining
#' singular values and corresponding singular vectors. 
#'
#' @param matr A numeric matrix to be pseudoinverted.
#'
#' @return
#' A numeric matrix giving the Moore-Penrose pseudoinverse of matr.
#' If all singular values are below the tolerance, a zero matrix of
#' appropriate dimensions is returned.
#'
#' @examples
#' \dontrun{
#' A <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
#' mp_inv(A)
#' }
#' 
#' @keywords internal
#' 
mp_inv <- function(matr) {
	svdm <- svd(matr)
	tolr <- max(dim(matr)) * .Machine$double.eps
    pstv <- (svdm$d > tolr * svdm$d[1])
    if (!any(pstv)) {
        return(matrix(0, ncol(matr), nrow(matr)))
    } else {
		temp <- svdm$v[, pstv, drop = FALSE] %*% diag(1/svdm$d[pstv], sum(pstv))
		return(tcrossprod(temp, svdm$u[, pstv, drop = FALSE]))
	}
}
