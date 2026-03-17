phzstat <- function(objt, krnl = "ba", band = "and", step = 2) {
	lgth <- objt$lgth
	yvls <- objt$yvls
	ylag <- objt$ylag
	deno <- objt$deno
	coef <- objt$coef
	rsds <- objt$rsds
	shrt <- objt$shrt
	cffs <- objt$cffs
	tsta <- objt$tsta
	#
    # Conventional, PP88
	#
    long <- lrvar(rsds, krnl, band)$longvar
    #
	zcff <- cffs - ((long - shrt)/2)/(deno/lgth^2)
	ztst <- sqrt(shrt/long) * tsta - ((long - shrt)/2)/sqrt(long * deno/lgth^2)
    #
	return(c(zcff = zcff, ztst = ztst))
}
