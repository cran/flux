## prepare function, making a list with a data.frame for
## each flux measurement (one chamber placement)
conz.prep <-
function(dat, columns = NULL, factors, nmes, min.cm = 3){
	## make the list for the by function
	sellist <- lapply(c(1:length(factors)), function(x) dat[,factors[x]])
	## prepare column vector in case columns = NULL
	if(is.null(columns)){ columns <- c(1:ncol(dat))}
	## separate into tables per chamber measurement
	conz.parts <- by(dat[,columns], sellist, function(x) x[])
	## in the prozess empty tables are produced when a combination
	## of factors in sellist has no data. the empty tables are
	## eliminated now
	flux.sel <- unlist(lapply(conz.parts, function(x) !is.null(x)))
	conz.parts <- conz.parts[flux.sel]
	## tables with less then a specified number of concentration
	## measurements (via min.cm) are eliminated as well
	flux.sel <- sapply(conz.parts, function(x) nrow(x)>=min.cm)
	conz.parts <- conz.parts[flux.sel]
	## do the naming (for easy access to the data)
	nams <- data.frame(t(sapply(conz.parts, function(x) sapply(x[][1,nmes], as.character))))
	nams$all <- apply(nams, 1, function(x) paste(sapply(x, as.character), sep=".", collapse="."))
	names(conz.parts) <- nams$all
	conz.parts <- list(tables = conz.parts, nmes = nams)
	return(conz.parts)
} 