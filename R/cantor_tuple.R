#' Cantor Tuple Function
#'
#' Extends the Cantor pairing function to tuples of arbitrary length.
#'
#' This function reduces an n-element numeric vector to a single non-negative
#' integer using repeated Cantor pairing. For n = 1 the value is returned
#' unchanged; for n = 2 it behaves like `cantor_pair`; for larger n it applies
#' recursive pairing.
#'
#' @param xval A numeric vector of length 1 or more.
#'
#' @return A single numeric value representing the Cantor-encoded tuple.
#' 
#' @seealso cantor_pair
#'
#' @examples
#' cantor_tuple(c(2, 4, 6))
#' 
#' @keywords internal
#' 
cantor_tuple <- function(xval) {
	lgth <- length(xval)
	if (lgth == 1) {
		return(xval)
	} else if (lgth == 2) {
		return(cantor_pair(xval))
	} else {
		return(cantor_pair(c(Recall(xval[1:(lgth - 1)]), xval[lgth])))
	}
}
