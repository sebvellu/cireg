getpozfbqntl <- function(
	prob, dpow, xnum, krnl, bfrc, lgth = 1000, simu = 10000,
	tolr = .Machine$double.eps
) {
	btbl <- 2 * 1:50/100
	btbl <- btbl[which.min(abs(btbl - bfrc))]
	tema <- try(lookupvewa23qntl(
		prob = prob,
		dpow = dpow,
		xnum = xnum,
		stat = "zcff",
		krnl = krnl,
		bfrc = btbl,
		tolr = tolr
	), TRUE)
	if (!inherits(tema, "try-error")) {
		if (!any(is.na(tema))) {
			temb <- lookupvewa23qntl(
				prob = prob,
				dpow = dpow,
				xnum = xnum,
				stat = "ztst",
				krnl = krnl,
				bfrc = btbl,
				tolr = tolr
			)
			return(drop(cbind(tema, temb, deparse.level = 0)))
		}
	}
	return(phzfbquant(prob, dpow, xnum, NULL, krnl, bfrc, lgth, simu))
}
