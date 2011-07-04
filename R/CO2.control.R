CO2.control <-
function(CO2.columns = c("CO2ppm", "CO2Code"), trans.out = TRUE, range.lim = 30, min.allow = 3){
	res <- list(columns = CO2.columns, trans.out = trans.out, range.lim = range.lim, min.allow = min.allow)
	return(res)
}