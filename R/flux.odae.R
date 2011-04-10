flux.odae <-
function(dat, min.allowed = 3, max.nrmse = 0.1, start = 1){
	ghg <- names(dat)[1]
	## const stays for constituent
	names(dat)[1:5] <- c("const", "time", "gc.qual", "volume", "temp")
	## transform data matrix according to the start flag
	## (which states, which zero measurement should be taken)
	dat.zero <- dat[dat$time == 0,]
	dat.runs <- dat[dat$time != 0,]
	dat <- rbind(dat.zero[start,], dat.runs)
	## getting all possible combinations of at least three x
	vers <- unlist(lapply(c(min.allowed:nrow(dat)), function(x) combn(c(1:nrow(dat)), x, simplify=FALSE)), recursive=FALSE)
	vers <- lapply(vers, "sort")
	## when 2 or more measurements were taken at a time the number
	## of min.allowed has to be extended to the number of unique
	## time steps. This is tested for and corrected in the following
	sel <- unlist(lapply(vers, function(x) length(unique(dat[x,2]))))
	vers <- vers[sel >= min.allowed]
	## determine number of concentration measurements per version
	vers.n <- sapply(vers, "length")
	## now we can calculate the regression models
	vers.lm <- lapply(vers, function(x) lm(const ~ time, data=dat[x,]))
	## and extract the statistic adj.r2
	#adj.r2 <- unlist(lapply(vers.lm, function(x) summary(x)$adj.r.squared))
	## as well as the statistic mrss (mean residual sum of squares)
	#mrss <- unlist(lapply(vers.lm, function(x) sqrt(sum(residuals(x)^2)/length(residuals(x)))))
	## and the statistic nrmse
	nrmse <- unlist(lapply(vers.lm, function(x) sqrt(sum(residuals(x)^2)/summary(x)$df[2])/diff(range(x$model[1], na.rm=TRUE))))
	## we rank them both
	#adj.r2.ranks <- order(adj.r2, decreasing = TRUE)
	#rss.ranks <- order(mrss)
	ranks <- order(vers.n, 1-nrmse, decreasing=TRUE)
	m2t <- ranks[nrmse[ranks] <= max.nrmse][1]
	if(is.na(m2t)){	
		m2t <- order(nrmse)[1]
	}
	lm4flux <- vers.lm[[m2t]]
	row.select <- vers[[m2t]]
	names(dat)[1] <- ghg
	res <- list(lm4flux = lm4flux, row.select = row.select, orig.dat=dat)
	return(res)
	}

