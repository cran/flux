inspect <-
function(x, what, sustain = NULL, retain = FALSE){
	z <- x$tables
	if(is.list(what)){
		if(!is.null(sustain)){names(what) <- paste(sustain, names(what), sep=".")}
		those <- match(names(what), names(z))
		for(i in c(1:length(what))){
			z[[those[i]]] <- z[[those[i]]][-what[[i]],]
		}
		if(retain){x$tables.orig <- x$tables}
		x$tables <- z
		return(x)
	}
	else{
		if(!is.null(sustain)){what <- paste(sustain, what, sep=".")}
		if(is.character(what)){what <- match(what, names(z))}
		return(z[what])
	}
}