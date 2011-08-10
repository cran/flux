"print.fluxes" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
	## prepare
	xt <- x$flux.table
	ab <- which(names(xt)=="all")
	nms.tab <- xt[,1:(ab-1)]
	CO2.tab <- data.frame(CO2.flags = paste(xt[,(ab+8)], xt[,(ab+6)], xt[,(ab+7)], ".", xt[,(ab+9)], xt[,(ab+10)], sep=""), CO2.flux = xt[,(ab+3)], CO2.us = format(paste(xt[,(ab+1)], xt[,(ab+2)], sep=" ")))
	CO2.tab[,2] <- round(CO2.tab[,2], 3)
	CH4.tab <- data.frame(CH4.flags = paste(xt[,(ab+18)], xt[,(ab+16)], xt[,(ab+17)], ".", xt[,(ab+19)], sep=""), CH4.flux = xt[,(ab+13)], CH4.us = format(paste(xt[,(ab+11)], xt[,(ab+12)], sep=" ")))
	CH4.tab[,2] <- round(CH4.tab[,2], 3)
	N2O.tab <- data.frame(N2O.flags = paste(xt[,(ab+27)], xt[,(ab+25)], xt[,(ab+26)], ".", xt[,(ab+28)], sep=""), N2O.flux = xt[,(ab+22)], N2O.us = format(paste(xt[,(ab+20)], xt[,(ab+21)], sep=" ")))
	N2O.tab[,2] <- round(N2O.tab[,2], 3)
	## rl stuff
	rla <- x$range.lim
	rl.out <- x$range.lim
	for(i in c(1:length(rla))){
		rl <- rla[[i]]
		if(length(rl)!=1){
			mean.rl <- round(mean(rl))
			range.rl <- round(range(rl))
			rl.statement <- paste(mean.rl, "on average, ranging from", range.rl[1], "to", range.rl[2])
		} else {rl.statement <- paste(rl, "(global)")}
		rl.out[[i]] <- rl.statement
	}
	cat("GHG flux rates and quality flags", "\n")
	cat("CO2 range limit:", rl.out[[1]], "\n")
	cat("CH4 range limit:", rl.out[[2]], "\n")
	cat("N2O range limit:", rl.out[[3]], "\n\n")
	## print to console
	print(cbind(CO2.tab, CH4.tab, N2O.tab))
	cat("\n\n")
	invisible(x)
}