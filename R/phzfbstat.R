phzfbstat <- function(objt, krnl = "ba", band = "and") {
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
    # Fixed-b asymptotic theory, VoWa13
    #
	cffm <- coef + (shrt/2)/(deno/lgth)
	rsdm <- yvls - cffm * ylag
	#
	lrvr <- lrvar(rsdm, krnl, band)
	fixb <- lrvr$bwdh/lgth
	long <- lrvr$longvar
	#
	zcfb <- cffs - ((long - shrt)/2)/(deno/lgth^2)
	ztsb <- sqrt(shrt/long) * tsta - ((long - shrt)/2)/sqrt(long * deno/lgth^2)
    #
	return(c(zcfb = zcfb, ztsb = ztsb, fixb = fixb))
}
