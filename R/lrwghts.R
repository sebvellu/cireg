#' Compute Lag Window Weights for Common Spectral Density Kernels
#'
#' This function generates lag window weights for several widely used
#' kernels in spectral density estimation. The available kernels are:
#' 
#' @param lgth Integer. Total length for which lag weights are generated.
#' 
#' @param bwdh Numeric. Bandwidth parameter controlling kernel scaling.
#' 
#' @param krnl Character string. Kernel type, one of:
#' 
#'   - `"ba"`: Bartlett
#'   - `"bo"`: Bohman
#'   - `"da"`: Daniell
#'   - `"pa"`: Parzen
#'   - `"qs"`: Quadratic Spectral
#'   - `"th"`: Tukey-Hanning
#'   - `"tr"`: Truncated
#' 
#' @note
#' Kernel weights are computed for lags from `1` up to the smaller of
#' `lgth - 1` or `ceiling(bwdh) - 1`, depending on the kernel.
#' 
#' @return
#' A numeric vector of kernel weights corresponding to the selected kernel.
#' 
#' @examples
#' \dontrun{
#' # Bartlett weights for length 20 with bandwidth 5
#' lrwghts(20, 5, "ba")
#'
#' # Quadratic spectral kernel example
#' lrwghts(30, 10, "qs")
#' }
#' 
#' @keywords internal
#' 
lrwghts <- function(lgth, bwdh, krnl = "ba") {
	if (krnl == "tr") { # Truncated/rectangular kernel
		return(rep(1, min(ceiling(bwdh) - 1, lgth - 1)))
	} else if (krnl == "ba") { # Bartlett kernel
		x <- seq_len(min(ceiling(bwdh) - 1, lgth - 1))/bwdh
		return(1 - x)
	} else if (krnl == "pa") { # Parzen kernel
		x <- seq_len(min(ceiling(bwdh) - 1, lgth - 1))/bwdh
		return(ifelse(x <= 1/2, 1 - 6 * x^2 + 6 * x^3, 2 * (1 - x)^3))
	} else if (krnl == "th") { # Tukey-Hanning kernel
		x <- seq_len(min(ceiling(bwdh) - 1, lgth - 1))/bwdh
		return((1 + cospi(x))/2)
	} else if (krnl == "qs") { # Quadratic spectral kernel
		x <- seq_len(lgth - 1)/bwdh
		x <- 6 * x/5
		return(3/(x * pi)^2 * (sinpi(x)/(pi * x) - cospi(x)))
	} else if (krnl == "bo") { # Bohman kernel
		x <- seq_len(min(ceiling(bwdh) - 1, lgth - 1))/bwdh
		return((1 - x) * cospi(x) + sinpi(x)/pi)
	} else { # if (krnl == "da") { # Daniell kernel
		x <- seq_len(lgth - 1)/bwdh
		return(sinpi(x)/(pi  * x))
	}
}
