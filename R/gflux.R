gflux <-
function(ct, c0, T, V, A = 1, M = 44, t = 1/60, p = 101325){
	R <- 8.314
	flux <- ((ct-c0) * V * M * p) / (t * R * (T + 273.15) * A)
	flux
	}

