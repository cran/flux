lips <- function(x, y, x.step = 1){
	sel <- (!is.na(x)) + (!is.na(y)) == 2
	x <- x[sel]
	y <- y[sel]
	x.diff <- x[-1] - x[-length(x)]
	if(class(x.diff)=="difftime"){
		x.diff <- difftime(x[-1], x[-length(x)], units="days")
		## check whether there is zero x.diff (must be treated special when x is a time)
		sel.x.diff0 <- x.diff==0
		if(sum(sel.x.diff0)>0){
			x[which(sel.x.diff0)+1] <- x[which(sel.x.diff0)+1] + x.step
			x.diff <- difftime(x[-1], x[-length(x)], units="days")
		}
		if(x.step!=1){
			x.diff <- as.numeric(x.diff*24*60*60/x.step)
		}
	}
	ns <- lapply(x.diff, function(z) seq(z)-1)
	if(is.factor(y)){
		y.out <- unlist(lapply(seq(length(ns)), function(x) rep(y[x], length(ns[[x]]))))
		y.out <- c(as.character(y.out), as.character(y.out[length(y.out)]))
	}
	else{
		y.diff <- y[-1] - y[-length(y)]
		slope <- y.diff/x.diff
		inc <- lapply(c(1:length(ns)), function(z) slope[[z]]*ns[[z]])
		y.out <- c(unlist(lapply(c(1:length(inc)), function(z) y[z]+inc[[z]])), y[length(y)])
	}
	x.out <- seq(min(x), max(x), x.step)
	res <- data.frame(x.out = x.out, y.out = y.out)
	return(res)
}