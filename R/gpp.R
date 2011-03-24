gpp <- function(NEE, PAR, PAR.Temp, Reco.m, allow.offset = FALSE){
	# predict Reco values using the Reco model provided
	Reco <- predict(Reco.m, newdata = data.frame(Temp = PAR.Temp))
	# calculate GPP
	GPP <- NEE - Reco
	# derive start values for the non-linear fitting
	# the start for alpha is taken from the linear regression
	# of GPP against PAR for PAR <= 500
	sel <- PAR <= 500
	coefs <- coef(lm(GPP ~ PAR, subset=sel))
	s.alpha <- coefs[2]
	# correct offset when wanted (default)
	offset <- 0
	if(!allow.offset){
		offset <- as.numeric(coefs[1])
		GPP <- GPP - offset
	}
	# the start value for GPmax is obtained by averaging the fluxes
	# at the 5 highest PAR values
	PAR.ord <- order(PAR)
	PAR.sel <- PAR.ord[(length(PAR.ord)-5) : length(PAR.ord)]
	s.GPmax <- mean(GPP[PAR.sel])
	# compile the start value list
	s.list <- list(GPmax = s.GPmax, alpha = s.alpha)
	# do the modeling
	mg <- nls(GPP ~ (GPmax * alpha * PAR)/(alpha * PAR + GPmax), start=s.list, model=TRUE)
	res <- list(mg = mg, mr = Reco.m, data = list(PAR.Temp = PAR.Temp, Reco = Reco, offset = offset))
	class(res) <- "gpp"
	return(res)
}