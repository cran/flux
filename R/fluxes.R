fluxes <-
function(x, var.par, co2ntrol = list(leak = TRUE, relay = FALSE), min.allowed = 3, max.nrmse = 0.1, nrmse.lim = 0.2, r2.qual = 0.8, range.lim = 30, out.unit = "auto", elementar = FALSE, hardflag = list(range = TRUE), p.values = TRUE){
	## extract the name vector from x
	nmes <- x$nmes
	## extract the data tables from x
	x <- x$tables
	## provide a leak.flag in case no co2ntrol is wanted
	leak.flag <- rep(2, length(x))
	## provide correct quality parameters for each gas
	if(length(max.nrmse)==1){max.nrmse <- list(CO2 = max.nrmse, CH4 = max.nrmse, N2O = max.nrmse)}
	if(length(nrmse.lim)==1){nrmse.lim <- list(CO2 = nrmse.lim, CH4 = nrmse.lim, N2O = nrmse.lim)}
	if(length(r2.qual)==1){r2.qual <- list(CO2 = r2.qual, CH4 = r2.qual, N2O = r2.qual)}
	if(length(range.lim)==1){range.lim <- list(CO2 = range.lim, CH4 = range.lim, N2O = range.lim)}
	if(length(out.unit)==1){out.unit <- list(CO2 = out.unit, CH4 = out.unit, N2O = out.unit)}
	## check about range.lim and attach calibration data to the data
		for(i in c(1:length(x))){
			x[[i]]$CO2.rl <- ifelse(length(range.lim$CO2)==1, range.lim$CO2, range.lim$CO2[i])
			x[[i]]$CH4.rl <- ifelse(length(range.lim$CH4)==1, range.lim$CH4, range.lim$CH4[i])
			x[[i]]$N2O.rl <- ifelse(length(range.lim$N2O)==1, range.lim$N2O, range.lim$N2O[i])
		}
	## run the outlier detection and elimination routines for the ghg
	## prepare selectors
	co <- grep("CO2", names(var.par))
	ch <- grep("CH4", names(var.par))
	no <- grep("N2O", names(var.par))
	cat(".",sep="")
	var.par.CO2 <- var.par
	names(var.par.CO2)[names(var.par.CO2)=="CO2"] <- "ghg"
	names(var.par.CO2)[names(var.par.CO2)=="CO2.gcq"] <- "gc.qual"
	CO2.pre <- lapply(x, function(x) flux.odae(x, var.par = var.par.CO2[-c(ch,no)], min.allowed = min.allowed, max.nrmse = max.nrmse$CO2, rl="CO2.rl"))
	cat(".",sep="")
	## do prep for CO2 and co2ntrol if wanted
	if(!is.null(co2ntrol)){
		CO2.range.check <- sapply(x, function(x) ifelse(diff(range(x[,var.par.CO2$ghg])) >= mean(x$CO2.rl), TRUE, FALSE))
		if(co2ntrol$leak){
			leak.flag <- sapply(CO2.pre, function(x) coef(x$lm4flux)[2]) <= 0
			leak.flag <- !as.logical(leak.flag*CO2.range.check)
		}
		if(co2ntrol$relay){
			x.sub <- lapply(which(CO2.range.check), function(y) x[[y]][CO2.pre[[y]]$row.select,])
			x[which(CO2.range.check)] <- x.sub
		}
	}
	cat(".",sep="")
	## do prep for CH4
	var.par.CH4 <- var.par
	names(var.par.CH4)[names(var.par.CH4)=="CH4"] <- "ghg"
	names(var.par.CH4)[names(var.par.CH4)=="CH4.gcq"] <- "gc.qual"
	names(var.par.CH4)[names(var.par.CH4)=="CH4.rl"] <- "rl"
	CH4.pre <- lapply(x, function(x) flux.odae(x, var.par = var.par.CH4[-c(co,no)], min.allowed = min.allowed, max.nrmse = max.nrmse$CH4, rl="CH4.rl"))
	cat(".",sep="")
	## do prep for N2O
	var.par.N2O <- var.par
	names(var.par.N2O)[names(var.par.N2O)=="N2O"] <- "ghg"
	names(var.par.N2O)[names(var.par.N2O)=="N2O.gcq"] <- "gc.qual"
	N2O.pre <- lapply(x, function(x) flux.odae(x, var.par = var.par.N2O[-c(co,ch)], min.allowed = min.allowed, max.nrmse = max.nrmse$N2O, rl="N2O.rl"))
	cat(".",sep="")
	## run the CO2 flux estimation via flux.conv and gflux
	CO2.res <- lapply(CO2.pre, function(x) flux.conv(x, ghg = "CO2", r2.qual = r2.qual$CO2, nrmse.lim = nrmse.lim$CO2, out.unit = out.unit$CO2, elementar = elementar, hardflag = hardflag))
	## run the CH4 flux estimation via flux.conv and gflux
	CH4.res <- lapply(CH4.pre, function(x) flux.conv(x, ghg = "CH4", r2.qual = r2.qual$CH4, nrmse.lim = nrmse.lim$CH4, out.unit = out.unit$CH4, elementar = elementar, hardflag = hardflag))
	## run the N2O flux estimation via flux.conv and gflux
	N2O.res <- lapply(N2O.pre, function(x) flux.conv(x, ghg = "N2O", r2.qual = r2.qual$N2O, nrmse.lim = nrmse.lim$N2O, out.unit = out.unit$N2O, elementar = elementar, hardflag = hardflag))
	## when pv values are to be reported… extract them
	if(p.values){
		CO2.pv <- sapply(CO2.res, function(x) coef(summary(x$fl.dat$lm4flux))[2,4])
		CO2.pv <- as.vector(symnum(CO2.pv, corr=FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," ")))
		CH4.pv <- sapply(CH4.res, function(x) coef(summary(x$fl.dat$lm4flux))[2,4])
		CH4.pv <- as.vector(symnum(CH4.pv, corr=FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," ")))
		N2O.pv <- sapply(N2O.res, function(x) coef(summary(x$fl.dat$lm4flux))[2,4])
		N2O.pv <- as.vector(symnum(N2O.pv, corr=FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," ")))
	}
	## make tables for the ghgases
	CO2.table <- t(sapply(CO2.res, function(x) unlist(x$fluss[2:8])))
	CO2.units <- sapply(CO2.res, function(x) x$unit)
	CO2.table <- data.frame(CO2.units, CO2.pv, CO2.table)
	CO2.table$leak.f <- leak.flag*1
	names(CO2.table)[-c(1:2)] <- paste("CO2.", names(CO2.table)[-c(1:2)], sep="")
	CH4.table <- t(sapply(CH4.res, function(x) unlist(x$fluss[2:8])))
	CH4.units <- sapply(CH4.res, function(x) x$unit)
	CH4.table <- data.frame(CH4.units, CH4.pv, CH4.table)
	names(CH4.table)[-c(1:2)] <- paste("CH4.", names(CH4.table)[-c(1:2)], sep="")
	N2O.table <- t(sapply(N2O.res, function(x) unlist(x$fluss[2:8])))
	N2O.units <- sapply(N2O.res, function(x) x$unit)
	N2O.table <- data.frame(N2O.units, N2O.pv, N2O.table)
	names(N2O.table)[-c(1:2)] <- paste("N2O.", names(N2O.table)[-c(1:2)], sep="")
	## handthrough data
	htd <- t(sapply(CO2.res, function(x) x$fl.dat$dat.out))
	## compile all in one big results table
	flux.table <- data.frame(nmes, CO2.table, CH4.table, N2O.table, htd[,-c(1:2)])
	## put the original results into a list
	flux.res <- list(CO2 = CO2.res, CH4 = CH4.res, N2O = N2O.res)
	## compile results for output
	res <- list(flux.res = flux.res, flux.table = flux.table, range.lim = range.lim)
	class(res) <- "fluxes"
	return(res)
}
