CO2.control <-
function(CO2.columns = c("CO2ppm", "CO2Code"), leak = TRUE, relay = FALSE, range.lim = 30, min.allow = 3, max.nrmse = 0.1){
	res <- list(columns = CO2.columns, leak = leak, relay = relay, range.lim = range.lim, min.allow = min.allow)
	return(res)
}