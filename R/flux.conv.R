flux.conv <-
function(fl.dat, ch.area, M = "CH4", range.lim = 30, r2.qual = 0.8, nrmse.lim = 0.2, nomba = NULL, in.unit = "ppm", out.unit = "mg"){
	lm4flux <- fl.dat$lm4flux
	dat <- fl.dat$orig.dat
	ct <- as.numeric(coef(lm4flux)[1] + coef(lm4flux)[2]*max(lm4flux$model[2]))
	c0 <- as.numeric(coef(lm4flux)[1])
	## set molar weight according to ghg
	m <- switch(M, CH4 = 16.04, N2O = 44.01, CO2 = 44.01)
	## average temperature in chamber during measurement
	T <- mean(dat$temp, na.rm=TRUE)
	V <- dat$volume[1]
	A <- ch.area
	t <- max(lm4flux$model[2])/60
	## calculating the flux with the Augustin equation
	## Details in gflux()!
	flux <- gflux(ct=ct, c0=c0, T=T, V=V, A=A, M=m, t=t)
	## according to the input unit the flux is transferred to the
	## common base g/m2*h
	flux <- switch(in.unit, ppm = flux/1e+6, ppb = flux/1e+9)
	## according to the output.unit the unit is changed
	## defaults to mg/m2*h
	flux <- switch(out.unit, ng = flux*1e+9, mug = flux*1e+6, mg = flux*1e+3, g = flux)
	## check the range limit and set flux to zero if below
	flux <- ifelse(diff(range(lm4flux$model[,1])) <= range.lim, 0, flux)
	## provide easy access to model r2.adj
	r2 <- summary(lm4flux)$r.squared
	## provide easy access to model normalized rss (normalized by averaging and squarerooting)
	#rss <- sqrt(sum(residuals(lm4flux)^2)/length(residuals(lm4flux)))
	## provide easy access to model nrmse
	nrmse <- sqrt(sum(residuals(lm4flux)^2)/summary(lm4flux)$df[2])/diff(range(lm4flux$model[1], na.rm=TRUE))
	## set the r2-quality flag
	r2.check <- ifelse(r2 >= r2.qual, TRUE, FALSE)
	r2.check <- ifelse(is.na(r2.check), FALSE, r2.check)
	## set the range-quality flag
	range.check <- ifelse(diff(range(lm4flux$model[,1])) >= range.lim, TRUE, FALSE)
	range.check <- ifelse(is.na(range.check), FALSE, range.check)
	## set the rss-quality flag
	#rss.check <- ifelse(rss <= rss.lim, TRUE, FALSE)
	#rss.check <- ifelse(is.na(rss.check), FALSE, rss.check)
	nrmse.check <- ifelse(nrmse <= nrmse.lim, TRUE, FALSE)
	nrmse.check <- ifelse(is.na(nrmse.check), FALSE, nrmse.check)
	## check nomba
	## ambient values from Mace Head Ireland and global average (CO2)
	## via http://cdiac.ornl.gov/pns/current_ghg.html
	## as of April 1st, 2011
	ambient <- switch(M, CH4 = 1870, N2O = 323, CO2 = 388.5)
	mba <- sum(dat[,1] <= ambient)
	if(is.null(nomba)){
		nomba <- ncol(dat)
	}
	nomba.check <- ifelse(mba <= nomba, TRUE, FALSE)
	nomba.check <- ifelse(is.na(nomba.check), FALSE, nomba.check)
	## prepare output
	fluss <- list(ghg = M, flux = flux, r2.check = r2.check, range.check = range.check, nrmse.check = nrmse.check, nomba.check = nomba.check, r2 = r2, nrmse = nrmse)
	## es fehlt noch die quality flag (nullfluss, r2check, rangecheck)
	res <- list(fluss = fluss, fl.dat = fl.dat, output.unit = out.unit)
	class(res) <- "flux"
	return(res)
	}

