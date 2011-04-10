flux <-
function(x, columns, M = "CH4", min.allowed = 3, max.nrmse = 0.1, nrmse.lim = 0.2, range.lim = 30, rl.backup = 30, r2.qual = 0.8, nomba = NULL, in.unit = "ppb", out.unit = "mg"){
	## extract the name vector from x
	nmes <- x$nmes
	## extract the data tables from x
	x <- x$tables
	## run the outlier detection and elimination routine (for details see flux.odae)
	flux.pre <- lapply(x, function(x) flux.odae(x[,columns], min.allowed = min.allowed, max.nrmse = max.nrmse))
	## run the flux calculation via flux.conv and gflux (consult these functions
	## for details)
	if(length(range.lim)>1){
		range.lim <- as.numeric(range.lim)
		range.lim[is.na(range.lim)] <- rl.backup
		flux.res <- lapply(c(1:length(flux.pre)), function(x) flux.conv(flux.pre[[x]], ch.area=flux.pre[[x]]$orig.dat[1,6], M = M, range.lim = range.lim[x], r2.qual = r2.qual, nrmse.lim = nrmse.lim, nomba = nomba, in.unit = in.unit, out.unit = out.unit))
	}
	else{
		flux.res <- lapply(flux.pre, function(x) flux.conv(x, ch.area=x$orig.dat[1,6], M = M, range.lim = range.lim, r2.qual = r2.qual, nrmse.lim = nrmse.lim, nomba = nomba, in.unit = in.unit, out.unit = out.unit))
	}
	## make table
	flux.table <- t(sapply(flux.res, function(x) unlist(x$fluss)))
	flux.table <- data.frame(nmes, flux.table)
	res <- list(flux.res = flux.res, flux.table = flux.table)
	class(res) <- "fluss"
	return(res)
}

