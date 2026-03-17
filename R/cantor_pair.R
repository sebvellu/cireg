#' Cantor Pairing Function
#'
#' Computes the Cantor pairing value for a 2-element numeric vector.
#'
#' The Cantor pairing function maps two non-negative integers (x, y)
#' to a single non-negative integer in a reversible way.
#'
#' @param xval A numeric vector of length 2 containing the values to pair.
#'
#' @return A single numeric value representing the Cantor pairing of `xval`.
#'
#' @examples
#' cantor_pair(c(3, 5))
#' 
#' @keywords internal
#' 
cantor_pair <- function(xval) {
	return((xval[1] + xval[2]) * (xval[1] + xval[2] + 1)/2 + xval[2])
}
