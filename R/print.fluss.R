"print.fluss" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
	tabl <- x$flux.table
	cat(as.character(tabl$ghg[1]), "flux rates", "\n\n")
	flags <- tabl[,c("r2.check", "range.check", "nrmse.check", "nomba.check", "leak.flag")]
	flags <- apply(flags, 2, function(x) as.logical(x)*1)
	flags <- data.frame(flags)
	flags$tight <- (flags$leak.flag == 0)*1
	names(flags) <- c("r2.f", "range.f", "nrmse.f", "nomba.f", "leak.f", "tight.f")
	num.vals <- apply(tabl[,c("flux", "r2", "nrmse")], 2, function(x) round(as.numeric(x), 3))
	ab <- which(names(tabl)=="all")
	new.table <- tabl[,1:(ab-1)]
	rownames(new.table) <- c(1:nrow(new.table))
	print(cbind(new.table, num.vals, flags[,c(1:4,6)]))
	cat("\n\n")
	invisible(x)
}