backsolid <- function(uppr) {# back solve identity
	return(backsolve(uppr, diag(1, ncol(uppr)))) #returns rinv
}