flux <-
function(x, columns, M = "CH4", co2ntrol = CO2.control(), min.allowed = 3, max.nrmse = 0.1, nrmse.lim = 0.2, range.lim = 30, rl.backup = 30, r2.qual = 0.8, nomba = NULL, in.unit = "ppb", out.unit = "mg"){
	## extract the name vector from x
	nmes <- x$nmes
	## extract the data tables from x
	x <- x$tables
	## provide an empty leak.flag vector (in case co2ntrol = NULL)
	leak.flag <- vector(length = length(x$tables))
	## check for necessity of co2ntrol and carry out...
	if(!is.null(co2ntrol)){
		columns.CO2 <- columns
		columns.CO2[c(1,3)] <- co2ntrol$columns
		## do CO2 range check on the original tables
		CO2.range.check <- sapply(x, function(x) ifelse(diff(range(x[,columns.CO2][,1])) >= co2ntrol$range.lim, TRUE, FALSE))
		## do flux.odae very much standard like on the CO2 data (have to be in x though)
		x.CO2 <- lapply(x, function(x) flux.odae(x[,columns.CO2], min.allowed = co2ntrol$min.allow, max.nrmse = max.nrmse))
		## calculate leak flag (CO2 is decreasing although an opaque chamber is used)
		leak.flag <- sapply(x.CO2, function(x) coef(x$lm4flux)[2]) <= 0
		leak.flag <- as.logical(leak.flag*CO2.range.check)
		## if methane and nitrous oxide concentration measurements shall be
		## skipped when the CO2 measurement at this time is an outlier
		if(co2ntrol$trans.out){
			x.sub <- lapply(which(CO2.range.check), function(y) x[[y]][x.CO2[[y]]$row.select,])
			x[which(CO2.range.check)] <- x.sub
		}
	}
	## run the outlier detection and elimination routine (for details see flux.odae)
	flux.pre <- lapply(x, function(x) flux.odae(x[,columns], min.allowed = min.allowed, max.nrmse = max.nrmse))
	## run the flux calculation via flux.conv and gflux (consult these functions
	## for details)
	if(length(range.lim)>1){
		range.lim <- as.numeric(range.lim)
		n.na <- sum(is.na(range.lim))
		if(n.na > 0){warning(paste(n.na, "of", nrow(nmes), "range limit values have been not available (NA)"))}
		range.lim[is.na(range.lim)] <- rl.backup
		flux.res <- lapply(c(1:length(flux.pre)), function(x) flux.conv(flux.pre[[x]], ch.area=flux.pre[[x]]$orig.dat[1,6], M = M, range.lim = range.lim[x], r2.qual = r2.qual, nrmse.lim = nrmse.lim, nomba = nomba, in.unit = in.unit, out.unit = out.unit))
		names(flux.res) <- nmes$all
	}
	else{
		flux.res <- lapply(flux.pre, function(x) flux.conv(x, ch.area=x$orig.dat[1,6], M = M, range.lim = range.lim, r2.qual = r2.qual, nrmse.lim = nrmse.lim, nomba = nomba, in.unit = in.unit, out.unit = out.unit))
	}
	## add leak flag to the flux.res fluss parts
	for(i in c(1:length(flux.res))){
		flux.res[[i]]$fluss$leak.flag <- leak.flag[i]
	}
	## make table
	flux.table <- t(sapply(flux.res, function(x) unlist(x$fluss)))
	flux.table <- data.frame(nmes, flux.table)
	res <- list(flux.res = flux.res, flux.table = flux.table)
	class(res) <- "fluss"
	return(res)
}

