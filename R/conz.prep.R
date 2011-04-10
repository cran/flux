## prepare function, making a list with a data.frame for
## each flux measurement (one chamber placement)
conz.prep <-
function(dat, columns, factors, nmes, min.cm = 3){
	## make the list for the by function
	sellist <- lapply(c(1:length(factors)), function(x) dat[,factors[x]])
	## separate into tables per chamber measurement
	conz.parts <- by(dat[,columns], sellist, function(x) x[])
	## in the prozess empty tables are produced when a comination
	## of factors in sellist has no data. the empty tables are
	## eliminated now
	flux.sel <- unlist(lapply(conz.parts, function(x) !is.null(x)))
	conz.parts <- conz.parts[flux.sel]
	## tables with less then a specified number of concentration
	## measurements (via min.cm) are eliminated as well
	flux.sel <- sapply(conz.parts, function(x) nrow(x)>=3)
	conz.parts <- conz.parts[flux.sel]
	## do the naming (for easy access to the data)
	names(conz.parts) <- sapply(conz.parts, function(x) paste(sapply(x[][1,nmes], as.character), sep=".", collapse="."))
	nams <- t(sapply(strsplit(names(conz.parts), "\\."), function(x) x[1:length(nmes)]))
	nams <- data.frame(nams)
	nams$all <- names(conz.parts)
	names(nams) <- c(nmes, "all")
	conz.parts <- list(tables = conz.parts, nmes = nams)
	return(conz.parts)
} 