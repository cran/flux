plot.reco <- function(x, ...){
	plot(R ~ Temp, data=x$model, ...)
	nd <- seq(min(x$model$Temp), max(x$model$Temp), length.out=100)
	lines(predict(x, newdata=data.frame(Temp=nd)) ~ nd)
}