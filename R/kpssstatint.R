kpssstatint <- function(rsds, long) {
	lgth <- length(rsds)
    return((sum(cumsum(rsds)^2)/long)/lgth^2)
}