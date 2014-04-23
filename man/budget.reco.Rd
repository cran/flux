\encoding{UTF-8}
\name{budget.reco}
\alias{budget.reco}
\alias{budget.gpp}
\title{
Predict fluxes from GPP and Reco models and prepare for summing them up to budgets.
}
\description{
The functions predict fluxes from GPP and \eqn{R_{eco}}{R_eco} models and prepare the data for summing them up to budgets including feeding in set.back positions and values (e.g. to truncate the output or to acknowledge the harvesting of biomass) and several adjustments and corrections.
}

\usage{
budget.reco(models, new.data, set.back = NULL, time.unit = "extract", 
adjust = TRUE, correct = list(thresh = "get", cvrm = TRUE, wndw = 12, 
intvl = 0.95), return.models = FALSE)

budget.gpp(models, new.data, set.back = NULL, time.unit = "extract", 
adjust = FALSE, correct = list(thresh = "get", cvrm = TRUE, wndw = 12, 
intvl = 0.95), PAR.correct = 100, return.models = FALSE)
}

\arguments{
  \item{models}{
List model objects of class \code{"breco"} or \code{"bgpp"}
}
  \item{new.data}{
Data from a climate station or other logging facility including timestamps. For \eqn{R_{eco}}{R_eco}  all temperature variables that have been used to construct the models are needed. For GPP PAR values are needed and the respective column must be named "PAR". The timestamp column (\code{POSIXt}) must be named "timestamp".
}
  \item{set.back}{
Data.frame with two columns, the first a \code{timestamp} and the second either a \code{value} of \eqn{R_{eco}}{R_eco} or GPP, to which fluxes should be set at the date and time in timestamp. Additionally, the second column can contain a value of \code{-999} or \code{-9999} for the start or end of the data series during budgeting.
}
  \item{time.unit}{
By default the time intervals between data points in \code{new.data} are extracted from the \code{timestamp} vector. However, by setting \code{time.unit} you can additionally give a number that represents the interval in seconds. The function uses this information to linearly interpolate values via \code{\link{lips}} when some are missing.
}
  \item{adjust}{
Logical. When \code{TRUE}, predicted fluxes are adjusted by fitting a linear model to the relationship of predicted to measured fluxes and extracting the slope. This  slope value is then used to adjust all predicted fluxes by calculating \code{flux} / \code{slope}. This is applied to the finished time series of fluxes. For \code{budget.reco} default is \code{TRUE} whereas for \code{budget.gpp} default is \code{FALSE}. See details.
}
  \item{correct}{
This triggers a further correction: Unreasonably high fluxes or spikes are identified and eliminated. This is applied to the finished time series of fluxes. The argument requires a named list with the entries \code{thresh}, \code{cvrm} (logical), \code{wndw}, and \code{intvl}. When no correction is wanted set \code{correct} = \code{NULL}. See details.
}
  \item{PAR.correct}{
Numeric representing the value of PAR (e.g. PPFD in micromole per squaremeter per second) below which GPP is assumed to be 0. This is applied to the finished time series of fluxes. All predicted GPP fluxes are corrected to 0 according to this value. See details.
}
  \item{return.models}{
Logical. Shall models be returned?
}
}

\details{
\bold{How does it work?} Both functions take a list of model objects (of class \code{"breco"} for \eqn{R_{eco}}{R_eco} models and of class \code{"bgpp"} for GPP) and predict fluxes via \code{\link[stats]{predict.nls}}. The required variables are taken from \code{new.data}. Fluxes are predicted for each model either for the whole time interval in \code{new.data} or when \code{set.back} is given with start and end configuration (values \code{-999} and \code{-9999}) for the defined time period. If the integration period is defined with \code{set.back}, the models that are temporally closest to the start and end times are duplicated and get the respective timestamps. With all models, fluxes are predicted forward from the date and time the models are hooked onto (argument \code{hook} in \code{\link{reco.bulk}} and \code{\link{gpp.bulk}}) up to the next model's date and time as well as backward down to the previous model's date and time. The resulting two flux vectors for the time period between two consecutive models (one based on the first, the other based on the second model) are linearly interpolated (see below). The same is done for the model errors. Errors are extracted from the model objects and assigned to each interpolated flux from one model to the next and to the previous. Then they are linearly interpolated (see  below).

Further, IDs are assigned to the weighted fluxes to identify the periods between two models. This is needed later for model error propagation. The fluxes in the time series between the first and the second model in the whole list get ID = 1, the fluxes in the time series between the second and the third model in the whole list get ID = 2, and so on.

\bold{How does interpolation work?} Linear interpolation of the two vectors (one based on the first model, the other based on the second model) is achieved by calculating weighted means of the integrated fluxes for each date and time with the weights the distances in time to the corresponding model timestamp. Thus, close to the one model timestamp the weighted mean of the two fluxes is almost entirely determined by that model, in the middle of the time series between two model dates both fluxes contribute equally to the mean and close to the other model timestamp the weighted mean of the two fluxes is almost entirely determined by that other model.

\bold{Why \code{adjust}?} It may happen that in sum the various models in a seasonal, annual or even bigger dat set tend to over- or underestimate the measured fluxes. To correct for this, the predicted fluxes can be adjusted as explained above. This typically leads to better overall modelling performance.

\bold{Why \code{correct}?} Especially with locally fitted \eqn{R_{eco}}{R_eco} models it happens - most often with winter data - that predicted fluxes are much higher than supported by the measurements because of the exponential element in fitting \eqn{R_{eco}}{R_eco} models. The model itself is fine but may have been fitted to temperature data spanning a relatively small range. If the temperatures in \code{new.data} are much higher (in relative terms), then unreasonable high fluxes may result. Such fluxes are identified and eliminated by \code{correct}. 

\bold{How is \code{correct}ed?} \code{thresh} specifies a maximum predicted flux allowed. When set to "get" (default) it is determined based on the data that were used to construct the models and represents the highest flux ever measured across all campaigns of the series. \code{cvrm} (logical) triggers whether further despiking should be done with \code{wndw} the width of the moving window in which the coefficient of variation (cv) is calculated and \code{intvl} the probability of the \code{\link[stats]{quantile}} against which cv should be tested. All fluxes with corresponding cv > \code{quantile(cv, intvl)} are eliminated.

\bold{Why \code{PAR.correct}?} Because the function not only predicts with all the provided models using \code{new.data} but also interpolates linearly between models, it happens that unreasonable GPP values are predicted reflecting photosynthesis under no or very low light conditions. These are just set to plant physiologically sensible 0. \code{PAR.correct}ion is done after inserting \code{set.back}(s) for cutting(s).

All data gaps resulting from corrections are then filled with linear interpolation from the values adjacent to the gap via \code{\link{lips}}.

\bold{How \code{set.back} is used to factor in cut dates or the like:}. When further \code{set.back} values are given in addition to the definition of the start and end dates, e.g. to acknowledge for biomass removal when predicting GPP fluxes, two things happen. First, the time series of fluxes resulting from all of the above is changed like this: At the date and time of a cut the flux  is set to the value given in \code{set.back} and then fluxes are linearly interpolated to the next proper model date by weighted means in the same manner as described above. Second, cut models are defined as linear models and integrated into the model list according to their timestamp. Thus, when run with \code{return.models = TRUE}, the updated model list can be used with \code{tbl8} to extract the relevant model parameters.
}

\value{
Both functions return a \code{data.frame} (called \code{tbl}) containing the predicted values, timestamps, etc. and optionally an object of class \code{"breco"} or \code{"bgpp"} (called \code{models}) containing the final models including the start and end models and any \code{set.back} models when \code{set.back}s were specified.

For \code{budget.reco} \code{tbl} has 4 columns
\item{reco.flux }{ Predicted and interpolated fluxes}
\item{reco.se }{ Model errors spawned across predicted and interpolated fluxes}
\item{reco.id }{ ID that identifies the model periods}
\item{timestamp }{ Timestamps}

For \code{budget.gpp} \code{tbl} has 4 columns
\item{gpp.flux }{ Predicted and interpolated fluxes}
\item{gpp.se }{ Model errors spawned across predicted and interpolated fluxes}
\item{gpp.id }{ ID that identifies the model periods}
\item{timestamp }{ Timestamps}

}

\references{
Beetz S, Liebersbach H, Glatzel S, Jurasinski G, Buczko U, Hoper H (2013) Effects of land use intensity on the full greenhouse gas balance in an Atlantic peat bog. Biogeosciences 10:1067-1082
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de},

with ideas by Sascha Beetz, \email{sascha.beetz@uni-rostock.de}
}

\seealso{
\code{\link{fluxx}}, \code{\link{reco}}, \code{\link{gpp}}, \code{\link{gpp2}}, \code{\link{reco.bulk}}, \code{\link{gpp.bulk}}, \code{\link{modjust}}
}

\examples{
### The examples are consecutive and are a suggestion
### how to run the whole process from bulk modelling
### over corrections and checks to full budgets including
### propagated model and interpolation error terms.

## The whole examples section is marked as
## not run because parts take longer than
## accepted by CRAN incoming checks.
## Remove first hash in each line to run them.

## Not run ##
## load data
#data(amd)
#data(amc)
## set global conversion factor
## All fluxes are in micromole * m-2 * s-1,
## thus they are transformed to g CO2-C * m-2 *1/2h
#uf <- 12*60*30/1000000
#
#### Reco ### (for details see reco.bulk)
## extract opaque (dark) chamber measurements 
#amr <- amd[amd$kind=="D",]
## fit reco models
#r.models <- reco.bulk(flux ~ t.air + t.soil2 + t.soil5 + 
#t.soil10 + timestamp, amr, amr$campaign, window=3, 
#remove.outliers=TRUE, method="arr", min.dp=2)
## adjust models
#r.models <- modjust(r.models, alpha=0.1, min.dp=3)
#
### prepare Reco budget (predict half hourly values for two years)
## define set.back with start and end dates only
## (typically you would take this from read in table)
#set.back <- data.frame(timestamp = c("2010-01-01 00:30", "2011-12-31 23:30"), 
#value = c(-999, -9999))
#set.back$timestamp <- strptime(set.back$timestamp, format="\%Y-\%m-\%d \%H:\%M", tz="GMT")
## run budget function with defaults
#r.bdgt <- budget.reco(r.models, amc, set.back)
## prepare for quality check (global model of predicted ~ measured)
#r.check <- checkm(r.bdgt, amr)
## have a look
#par(pty="s")
#lims <- range(r.check$reco.flux, r.check$flux)
#plot(reco.flux ~ flux, data=r.check, xlim=lims, ylim=lims)
#abline(coef=c(0,1), lty=3)
#mf1 <- lm(reco.flux ~ flux, r.check)
#abline(mf1)
#summary(mf1)
## calculate daily values (better for plotting)
#r.bdgt$day <- format(r.bdgt$timestamp, format="\%Y-\%m-\%d")
#r.daily <- data.frame(day = unique(r.bdgt$day))
#r.daily$day <- as.POSIXct(strptime(r.daily$day, format="\%Y-\%m-\%d", 
#tz="GMT"), tz="GMT")
## in addition to summing up per day, change unit
#r.daily$reco <- tapply(r.bdgt$reco.flux*uf, r.bdgt$day, sum)
## same for the model error terms
#r.daily$se <- tapply(r.bdgt$reco.se*uf, r.bdgt$day, sum)
#
#
#### GPP ### (for details see gpp.bulk)
#g.models <- gpp.bulk(flux ~ PAR + timestamp + kind, amd, amd$campaign, 
#method="Falge", min.dp=5)
#
### prepare GPP budget (predict half hourly values for two years)
## define set.back with start, end, and cut dates
## (typically you would take this from read in table)
#set.back <- data.frame(timestamp = c("2010-01-01 00:30", "2011-12-31 23:30", 
#"2010-07-22 12:00", "2010-09-03 12:00", "2010-10-13 12:00", "2011-06-29 12:00", 
#"2011-08-11 12:00", "2011-10-21 12:00"), value = c(-999, -9999, rep(-0.0001, 6)))
#set.back$timestamp <- strptime(set.back$timestamp, format="\%Y-\%m-\%d \%H:\%M", tz="GMT")
## have a look at the resulting data.frame
#set.back
## run budget function with correct = NULL and return the models
#g.bdgt <- budget.gpp(g.models, amc, set.back, correct=NULL, return.models=TRUE)
## the cut models are also returned:
#tbl8(g.bdgt$models)
## extract the half hourly values to proceed
#g.bdgt <- g.bdgt$tbl
## make daily values (better for plotting)
#g.bdgt$day <- format(g.bdgt$timestamp, format="\%Y-\%m-\%d")
#g.daily <- data.frame(day = unique(g.bdgt$day))
#g.daily$day <- as.POSIXct(strptime(g.daily$day, format="\%Y-\%m-\%d", tz="GMT"), tz="GMT")
## in addition to summing up per day, change unit
#g.daily$gpp <- tapply(g.bdgt$gpp.flux*uf, g.bdgt$day, sum)
## same for the model error terms
#g.daily$se <- tapply(g.bdgt$gpp.se*uf, g.bdgt$day, sum)
#
#
#### Budgets ###
### doing the actual budgeting
## first bring Reco and GPP budget data together
## because of different handling data sets
## may be of different length, therefore use merge
#r.bdgt$ts <- as.character(r.bdgt$timestamp)
#g.bdgt$ts <- as.character(g.bdgt$timestamp)
#bdgt <- merge(r.bdgt, g.bdgt[, -c(ncol(g.bdgt)-2, ncol(g.bdgt)-1)], 
#by.x="ts", by.y="ts", all.x=TRUE)
## calculate NEE
#bdgt$nee.flux <- bdgt$reco.flux + bdgt$gpp.flux
## error propagation 
#bdgt$nee.se <- sqrt(bdgt$reco.se^2 + bdgt$gpp.se^2)/sqrt(2)
## define unique id that spans across reco and gpp ids
#bdgt$nee.id <- paste(bdgt$reco.id, bdgt$gpp.id, sep=".")
## do budgets of fluxes (sum and use global conversion factor, see above) and error terms
## the model errors are summed up per model id and resulting
## sums are combined following error propagation
#with(bdgt, {c(
#  reco = sum(reco.flux*uf, na.rm=TRUE), 
#	reco.me = sqrt(sum(tapply(reco.se, reco.id, sum)^2))*uf, 
#	gpp = sum(gpp.flux*uf, na.rm=TRUE), 
#	gpp.me = sqrt(sum(tapply(gpp.se, gpp.id, sum)^2))*uf, 
#	nee = sum(nee.flux*uf, na.rm=TRUE), 
#	nee.me = sqrt(sum(tapply(nee.se, nee.id, sum, na.rm=TRUE)^2))*uf
#)})
#
### annual budget incl. interpolation error
#set.back <- data.frame(timestamp = c("2010-01-01 00:30", 
#"2010-12-31 23:30"), value = c(-999, -9999))
#set.back$timestamp <- strptime(set.back$timestamp, format="\%Y-\%m-\%d \%H:\%M", tz="GMT")
### reco.ie
## redoing budget with annual bounds
#r.bdgt.2010 <- budget.reco(r.models, amc, set.back, return.models=TRUE)
## run budget.ie with lo = 3 and it = 10 (default of 100 is advisable but slow)
#r.ie2010 <- budget.ie(r.bdgt.2010, lo=3, it=10) 
## summary statistic and factor in the uf
#r.ie2010 <- sd(r.ie2010*uf) 
## gpp.ie
## redoing budget with annual bounds
#g.bdgt.2010 <- budget.gpp(g.models, amc, set.back, correct=NULL, return.models=TRUE)
#g.ie2010 <- budget.ie(g.bdgt.2010, lo=3, it=10) 
#g.ie2010 <- sd(g.ie2010*uf) 
#
## do the actual budgeting
#tmp <- bdgt[(bdgt$timestamp >= set.back$timestamp[1]) & 
#(bdgt$timestamp <= set.back$timestamp[2]),]
#with(tmp, {c(
#  reco = sum(reco.flux*uf, na.rm=TRUE), 
#	reco.me = sqrt(sum(tapply(reco.se, reco.id, sum)^2))*uf,
#	reco.ie = r.ie2010,
#	gpp = sum(gpp.flux*uf, na.rm=TRUE), 
#	gpp.me = sqrt(sum(tapply(gpp.se, gpp.id, sum)^2))*uf,
#	gpp.ie = g.ie2010, 
#	nee = sum(nee.flux*uf, na.rm=TRUE), 
#	nee.me = sqrt(sum(tapply(nee.se, nee.id, sum, na.rm=TRUE)^2))*uf,
#	nee.ie = sqrt(r.ie2010^2 + g.ie2010^2)
#)})
#
## End not run ##
}

\keyword{ manip }
\keyword{ univar }
