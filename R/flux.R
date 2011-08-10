flux <-
function(x, var.par, ghg = "CH4", co2ntrol = CO2.control(), min.allowed = 3, max.nrmse = 0.1, nrmse.lim = 0.2, r2.qual = 0.8, range.lim = 30, out.unit = "auto", elementar = FALSE, hardflag = list(range = TRUE), p.values=TRUE){
	## get call for later reference
	cl <- match.call()
	## extract the name vector from x
	nmes <- x$nmes
	## extract the data tables from x
	x <- x$tables
	## provide an empty leak.flag vector (in case co2ntrol = NULL)
	leak.flag <- rep(2, length(x))
	## check about range.lim and attach calibration data to the data
	## and provide for handing through the correct range.lim to the results
	if(!is.null(range.lim)){
		for(i in c(1:length(x))){
			x[[i]]$rl <- ifelse(length(range.lim)==1, range.lim, range.lim[i])
		}
		if(is.character(range.lim)){
			range.lim <- sapply(x, function(x) x[1,range.lim])
		}
	}
	else{
		if(exists("rl", x[[1]])){
			range.lim <- sapply(x, function(x) mean(x$rl))
		} else {range.lim <- 0}
	}
	## check if co2ntrol is wanted and carry out when it is…
	if(!is.null(co2ntrol)){
		## check about range.lim and attach calibration data to the data
		for(i in c(1:length(x))){
			x[[i]]$CO2.rl <- ifelse(length(co2ntrol$range.lim)==1, co2ntrol$range.lim, co2ntrol$range.lim[i])
		}
		var.par.CO2 <- var.par
		var.par.CO2[c("ghg", "gc.qual")] <- co2ntrol$columns
		## do CO2 range check on the original tables
		CO2.range.check <- sapply(x, function(x) ifelse(diff(range(x[,var.par.CO2$ghg])) >= mean(x$CO2.rl), TRUE, FALSE))
		## do flux.odae very much standard like on the CO2 data (have to be in x though)
		x.CO2 <- lapply(x, function(x) flux.odae(x, var.par = var.par.CO2, min.allowed = co2ntrol$min.allow, max.nrmse = co2ntrol$max.nrmse, rl = "CO2.rl"))
		## calculate leak flag (CO2 is decreasing although an opaque chamber is used)
		if(co2ntrol$leak){
			leak.flag <- sapply(x.CO2, function(x) coef(x$lm4flux)[2]) <= 0
			leak.flag <- !as.logical(leak.flag*CO2.range.check)
		}
		## if methane and nitrous oxide concentration measurements shall be
		## skipped when the CO2 measurement at this time is an outlier
		if(co2ntrol$relay){
			x.sub <- lapply(which(CO2.range.check), function(y) x[[y]][x.CO2[[y]]$row.select,])
			x[which(CO2.range.check)] <- x.sub
		}
	}
	## run the outlier detection and elimination routine (for details see flux.odae)
	flux.pre <- lapply(x, function(x) flux.odae(x, var.par = var.par, min.allowed = min.allowed, max.nrmse = max.nrmse))
	## run the flux calculation via flux.conv and gflux (consult these functions
	## for details)
	flux.res <- lapply(flux.pre, function(x) flux.conv(x, ghg = ghg, r2.qual = r2.qual, nrmse.lim = nrmse.lim, out.unit = out.unit, elementar = elementar, hardflag = hardflag))
	## add leak flag to the flux.res fluss parts
	for(i in c(1:length(flux.res))){
		flux.res[[i]]$fluss$leak.f <- leak.flag[i]
	}
	## make table
	flux.table <- t(sapply(flux.res, function(x) unlist(x$fluss[1:8])))
	units <- sapply(flux.res, function(x) x$unit)
	flux.table <- data.frame(nmes, unit = units, flux.table, leak.f = leak.flag)
	if(p.values){
		pv <- sapply(flux.res, function(x) coef(summary(x$fl.dat$lm4flux))[2,4])
		flux.table$pv <- pv
	}
	## hand through data
	htd <- t(sapply(flux.res, function(x) x$fl.dat$dat.out))
	flux.table <- cbind(flux.table, htd[,-c(1:2)])
	res <- list(call = cl, flux.res = flux.res, flux.table = flux.table, range.lim = range.lim)
	class(res) <- "fluss"
	return(res)
}

