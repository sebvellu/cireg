bwclc <- function(tsrs, krnl = "ba", band = "and", wght = NULL) {
	if (identical(band, "and")) { # calculate Andrews 1991 bandwidth:
		tsrs <- as.matrix(tsrs)
		lgth <- nrow(tsrs)
		nums <- ncol(tsrs)
		if (is.null(wght)) {
			wght <- rep(1, nums)
		}
		estm <- matrix(nrow = 2, ncol = nums) 
		for (i in 1:nums) {
			yvls <- tsrs[2:lgth, i]
			xvls <- tsrs[1:(lgth - 1), i]
			estm[1, i] <- sum(yvls * xvls)/sum(xvls^2)
			estm[2, i] <- mean((yvls - estm[1, i] * xvls)^2)
		}
		# Andrews 1991, p. 834 (6.2) [and footnote 5] and p. 835 (6.4)
		if (krnl == "ba") {
			band <- (1 - estm[1, ])^6 * (1 + estm[1, ])^2
			band <- wght * 4 * estm[1, ]^2 * estm[2, ]^2/band
			band <- sum(band)/sum(wght * estm[2, ]^2/(1 - estm[1, ])^4)
			return(1.1447 * (band * lgth)^(1/3))
		} else {
			band <- (1 - estm[1, ])^8
			band <- wght * 4 * estm[1, ]^2 * estm[2, ]^2/band
			band <- sum(band)/sum(wght * estm[2, ]^2/(1 - estm[1, ])^4)
			if (krnl == "tr") {
				return(0.6611 * (band * lgth)^(1/5))
			} else if (krnl == "pa") {
				return(2.6614 * (band * lgth)^(1/5))
			} else if (krnl == "th") {
				return(1.7462 * (band * lgth)^(1/5))
			} else { #if (krnl == "qs") {
				return(1.3221 * (band * lgth)^(1/5))
			}
		}
	} else if (identical(band, "nw")) { # calculate Newey and West 1994
		# Newey and West 1994, p. 640--641, Table I, II (non-prewhitened):
		tsrs <- as.matrix(tsrs)
		lgth <- nrow(tsrs)
		if (is.null(wght)) {
			wght <- rep(1, ncol(tsrs))
		}
		if (krnl == "ba") {
			keys <- c(1, 1.1447, floor(4 * (lgth/100)^(2/9)))
		} else if (krnl == "pa") {
			keys <- c(2, 2.6614, floor(4 * (lgth/100)^(4/25)))
		} else { #if (krnl == "qs") {
			keys <- c(2, 1.3221, floor(4 * (lgth/100)^(2/25)))
		}
		tsrs <- drop(tsrs %*% wght)/sqrt(lgth)
		asum <- c()
		for (i in 0:keys[3]) {
			asum[i + 1] <- sum(tsrs[(i + 1):lgth] * tsrs[1:(lgth - i)])
		}
		dnom <- asum[1] + 2 * sum(asum[-1])
		numr <- 2 * sum((1:keys[3])^keys[1] * asum[-1])
		band <- keys[2] * ((numr/dnom)^2)^(1/(2 * keys[1] + 1))
		# 
		return(band * lgth^(1/(2 * keys[1] + 1)))
	} else if (identical(band, "nwt")) {
		return(floor(4 * (T/100)^(2/9)))
	} else {
		return(band)
	}
}
