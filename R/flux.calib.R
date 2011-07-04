flux.calib <-
function(dat, columns, calib, format = "%d.%m.%Y", calib.max = 12000, range.ext = 200, window = 7, calib.gas.defaults = c(300, 1000, 2000, 6000, 10000)){
	# defining the function which does the work
	flux.cal <-
	function(conz.dat, calib, format = "%d.%m.%Y", calib.max=12000, range.ext=200, window=7, calib.gas.defaults=c(300, 1000, 2000, 6000, 10000)){
		# extract date
		dts.dat <- as.Date(strptime(conz.dat[,1], format = format))
		m.date <- dts.dat[1]
		dts.cal <- as.Date(strptime(calib[,1], format = format))
		# extract conzentration data and calculate range plus padding (giving by 
		# the argument range.ext in the same units as the concentration measurements)
		conz.range <- range(conz.dat[,2], na.rm=TRUE) + c(-range.ext, range.ext)
		# skip calibration gas measurements with 
		calib <- calib[calib[,2] < calib.max,]
		# extract calibration gas measurements according to the date of the measurement
		# of the ghg and a window width window (days) around it
		calib <- calib[(dts.cal >= (m.date-window)) & (dts.cal <= (m.date+window)), 2]
		# extract only those calibration gases that fall into the conzentration range
		# its a bit more than the pure range (remember, above we added a buffer)
		calib.sub <- calib[(calib > conz.range[1]) & (calib < conz.range[2])]
		# calculate range limits (standard deviation of the calibration gas 
		# measurements) per calibration gas  
		range.lims <- tapply(calib.sub, cut(calib.sub, c(0, runmean(calib.gas.defaults, 2, endrule="trim")), labels=FALSE), sd)
		# range.means <- tapply(calib.sub, cut(calib.sub, c(0,runmean(calib.gas.defaults, 2, endrule="trim")), labels=FALSE), mean)
		# calculate average range limits across all included calibration gases
		range.lim <- mean(range.lims)
		return(range.lim)
	}

	# actually do the work
	# extract the needed columns from calib
	calib <- calib[,columns]
	ghg.lim <- sapply(dat$tables, function(x) flux.cal(x[,columns], calib[,columns], format = format, calib.max = calib.max, range.ext = range.ext, window = window, calib.gas.defaults = calib.gas.defaults))
	return(ghg.lim)
}