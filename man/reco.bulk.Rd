\encoding{UTF-8}
\name{reco.bulk}
\alias{reco.bulk}
\alias{gpp.bulk}
\title{
Bulk fitting of Reco and GPP models
}
\description{
The function allows for bulk fitting of \eqn{R_{eco}}{R_eco} and GPP models with the respective functions \code{\link{reco}} and \code{\link{gpp}}. This is often appropriate because data are gathered over a season, a year or longer...
}

\usage{
reco.bulk(formula, data, INDEX, window = 1, hook = "mean", remove.outliers = FALSE, 
fall.back = TRUE, ...)

gpp.bulk(formula, data, INDEX, window = 1, hook = "mean", oot.id = c("D", "T"), 
min.dp = 5, Reco.m = NULL, ts.Reco = NULL, fall.back = TRUE, ...)
}

\arguments{
  \item{formula}{
An object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the terms that are used in bulk \eqn{R_{eco}}{R_eco} and GPP model fitting. Choices of terms are more restricted than typically (see details). For instance, a timestamp always has to be provided. Also, temperature variables are required for \code{gpp.bulk} if \eqn{R_{eco}}{R_eco} values are predicted from models.
}
  \item{data}{
A data frame (or an object that can be coerced to that class by \code{as.data.frame}) containing at least all the 'model' terms specified in \code{formula}.
}
  \item{INDEX}{
A vector of length \code{nrow(data)} that is used to extract and compile data, for instance according to measurement campaign in the field. Internally \code{\link[base]{split}} is used with \code{f = INDEX} to create a list of \code{data.frame}s of which each contains all flux measurements for one model.
}
  \item{window}{
Both functions can fit the respective models across a moving window of adjacent \code{INDEX} values. Not advisable for GPP while \eqn{R_{eco}}{R_eco} modelling can really profit because more data points often lead to better models.
}
  \item{hook}{
Character string specifiying the kind of summary statistics used to fix a date and time to which  the fitted model shall refer. Up to now this is simply achieved by doing one of these summary statistics on the timestamp: \code{\link[base]{mean}}, \code{\link[base]{min}}, \code{\link[base]{max}} or \code{\link[stats]{median}}.
}
  \item{remove.outliers}{
Logical. If \code{TRUE} the function searches for outliers in the data points of the \eqn{R_{eco}}{R_eco} models and eliminates them. Per model the \code{\link[grDevices]{boxplot.stats}} of the residuals are obtained and if outliers are present they are eliminated and the model is fitted again. This is done twice. If the function fails in fitting the model to the new data set it falls back to the original data.
}
  \item{oot.id}{
Vector of length 2 that specifies which of the flux values derive from opaque (first value, i.e. \eqn{R_{eco}}{R_eco} measurements) and which derive from transparent (second value, i.e. NEE measurements) chamber measurements when \code{data} contains both. This is one of several approaches to GPP modeling here. See details.
}
  \item{min.dp}{
Numeric. Specifies the minimum number of data points that are accepted per model. Defaults to 5 which is already quite a small number.
}
  \item{Reco.m}{
Either an object of class "\code{reco}" resulting from \code{\link{reco}} with one \eqn{R_{eco}}{R_eco} model or an object of class "\code{breco}" resulting from \code{\link{reco.bulk}} with several \eqn{R_{eco}}{R_eco} models or a vector with estimated half hourly or hourly (or whatever interval you have) \eqn{R_{eco}}{R_eco} values. In case of the latter \code{ts.Reco} has to be specified as well because it is also used as a switch between internal \eqn{R_{eco}}{R_eco} modeling and assigning existing \eqn{R_{eco}}{R_eco} values. See details.
}
  \item{ts.Reco}{
POSIXlt or POSIXct vector with timestamps of the fluxes in \code{Reco.m}. Further, the default (\code{ts.Reco} = NULL) lets the function expect model object(s) in \code{Reco.m}.
}
  \item{fall.back}{
Logical. When TRUE the function falls back to linear mean models when the non linear approach did not work out (for \code{reco.bulk}: the slope of the linear relationship between Reco and temperature is < 0; for \code{gpp.bulk}: either no model could be fit or the starting slope parameter \code{alpha} is > 0). To do so a virtual data set is created with 50 random \code{GPP} values that have the same mean and sd as the original data and with a sequence of 50 \code{PAR} values spanning from 0 to 2000. A linear model is fit to these data with \code{lm(GPP ~ PAR)}.
}
  \item{\dots}{
Further arguments passed to \code{\link{reco}} or \code{\link{gpp}} e.g., the method for fitting the model when not using the respective defaults.
}
}

\details{
Models are - comparable to regression models - specified symbolically. Accordingly, the basic form is \code{response ~ terms} with \code{response} always referring to \eqn{CO_2}{CO2} exchange rates. For \code{terms} requirements differ between the two methods. In contrast to other \code{formula}e the \code{response} and all \code{terms} have to be in \code{data}.

\code{reco.bulk} expects a \code{formula} of the form \code{Reco} ~ \code{T1} +  \code{...} + \code{timestamp} with \code{Reco} referring to \eqn{CO_2}{CO2} fluxes estimated based on opaque chamber measurements (for instance with \code{\link{flux}}), \code{T1} referring to temperature readings relevant for \code{Reco} (e.g. air temperature) and taken during the corresponding chamber measurements. The \code{...} symbolizes that several more temperature readings can be specified if available (e.g. temperature in soil at 2cm), as many as you want. When more than one temperature is specified models are fit for each temperature and the best one is determined via \code{\link[stats]{AIC}} and reported together with the name of the corresponding temperature variable. Finally, \code{timestamp} is referring to the \code{POSIXt} timestamps that represent the dates and times of the corresponding measurements. \code{timestamp} always has to be specified as the last term of the \code{formula}. Models are fit using \code{\link{reco}}.

\code{gpp.bulk} expects a \code{formula} of the form \code{NEE} ~ \code{PAR} +  \code{timestamp} + \code{...} with \code{NEE} referring to \eqn{CO_2}{CO2} fluxes estimated based on transparent chamber measurements (for instance with \code{\link{flux}}), \code{PAR} referring to readings of the photosynthetically active radiation relevant for \code{NEE} and taken during the corresponding chamber measurements. The \code{...} symbolizes that several more terms can or have to be specified. This depends on the approach to the \eqn{R_{eco}}{R_eco} part of the GPP modeling (see \code{\link{gpp}}).

Approaches to estimate GPP values from measured NEE data using corresponding \eqn{R_{eco}}{R_eco} values:

\bold{Approach 1:} Extract corresponding \eqn{R_{eco}}{R_eco} fluxes from the provided data that are assigned to corresponding NEE values via their timestamp: For this approach \code{data} has to contain both NEE and \eqn{R_{eco}}{R_eco} fluxes and the model formula is specified as \code{NEE} ~ \code{PAR} +  \code{timestamp} + \code{oot} with the latter referring to a variable that indicates whether the respective fluxes were measured as \code{NEE} (transparent chamber) or \code{Reco} (opaque chamber or low PAR). In addition \code{oot.id} may have to be changed accordingly. \code{\link{gpp2}} is used for fitting the models. 

\bold{Approach 2:} Provide measured \eqn{R_{eco}}{R_eco} fluxes that are assigned to corresponding NEE values via their timestamp: To do this set \code{ts.Reco} != \code{NULL} and \code{Reco.m} a vector of \eqn{R_{eco}}{R_eco} fluxes and specifiy model with: \code{NEE} ~ \code{PAR} +  \code{timestamp}. \code{\link{gpp}} is used for fitting the models.

\bold{Approach 3:} Provide one \eqn{R_{eco}}{R_eco} model to predict \eqn{R_{eco}}{R_eco} fluxes at the time of the NEE measurements using the same temperature variable that was used to construct the model (with \code{\link{reco}}). Specify model with: \code{NEE} ~ \code{PAR} +  \code{timestamp} + \code{temperature}. \code{\link{gpp}} is used for fitting the models.

\bold{Approach 4:} Provide several \eqn{R_{eco}}{R_eco} models to predict \eqn{R_{eco}}{R_eco} fluxes at the time of the NEE measurements using the same temperature variables that were used to construct the models (with \code{\link{reco.bulk}}). The corresponding models are assigned to the NEE data via the timestamps that they carry. Specify model with: \code{NEE} ~ \code{PAR} +  \code{timestamp} + \code{temperature1} + \code{temperature2} + \code{temperature3} + \code{...}. All temperatures that may have been used for fitting the \eqn{R_{eco}}{R_eco} models (see above) should be given. \code{\link{gpp}} is used for fitting the models.

\code{remove.outliers} may result in better \eqn{R_{eco}}{R_eco} models. One should be careful with this and watch out for cases in which too many data points are eliminated. The function returns the number of skipped outliers per model to do just that.

If \code{fall.back} = \code{TRUE} no failed model fits are reported. That's quite useful when further bulk methods like \code{\link{budget.reco}} or \code{\link{budget.gpp}} are used to get annual or seasonal budgets.
}

\value{
Both functions return complex list structures with models. 

Output of \code{reco.bulk}:
Object of class "\code{breco}", a list with \code{length(unique(INDEX))} elements, each containing 3 elements:
\item{ts }{ Timestamp of the model.}
\item{mod }{ Has itself two elements. The first contains the model object as returned by \code{reco} and is named according to the method used. The second, \code{n.out}, is optional (only reported when \code{remove.outliers = TRUE} and there were indeed outliers identified and skipped) and gives the number of omitted data points.}
\item{which.Temp }{ Character string that identifies the temperature variable that was finally used for constructing the best model}

Output of \code{gpp.bulk}:
Object of class "\code{bgpp}", a list with \code{length(unique(INDEX))} elements each containing itself 2 entries:
\item{ts }{ Timestamp of the model}
\item{mod }{ Either an object of class "\code{gpp}" or of class "\code{gpp2}" depending on the approach used. Approaches 1 and 2 return "\code{gpp2}" objects, Approaches 3 and 4 return "\code{gpp}" objects. See \code{\link{gpp}} and \code{\link{gpp2}} for details.}
}

\references{
Beetz S, Liebersbach H, Glatzel S, Jurasinski G, Buczko U, Hoper H (2013) Effects of land use intensity on the full greenhouse gas balance in an Atlantic peat bog. Biogeosciences 10:1067-1082
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de},

with suggestions by Sascha Beetz, \email{sascha.beetz@uni-rostock.de}
}

\seealso{
\code{\link{reco}}, \code{\link{gpp}}, \code{\link{gpp2}}, \code{\link{fluxx}}, \code{\link{modjust}}
}

\examples{
## Whole example is consecutive and largely marked as
## not run because parts take longer than
## accepted by CRAN incoming checks.
## Remove first hash in each line to run them.
data(amd)
data(amc)

### Reco ###
## do reco models with 3 campaign wide window and 
## outlier removal (outliers according to models)
# first extract opaque (dark) chamber measurements 
amr <- amd[amd$kind=="D",]

## Nor run ##
## do bulk fitting of reco models (all specified temperatures 
## are tested and the best model (per campaign) is finally stored)
#r.models <- reco.bulk(flux ~ t.air + t.soil2 + t.soil5 + 
#t.soil10 + timestamp, amr, amr$campaign, window=3, 
#remove.outliers=TRUE, method="arr", min.dp=2)
#
## adjust models (BEWARE: stupid models with t1 >= 20 are skipped 
## within the function, this can be changed)
#r.models <- modjust(r.models, alpha=0.1, min.dp=3)
#
## make data.frame (table) for overview of model parameters
## the temperature with which the best model could be fit is reported
## this information also resides in the model objects in r.models
#tbl8(r.models)
#
#### GPP ###
### fit GPP models using method = Falge and min.dp = 5
### and take opaque (dark, i.e. reco) measurements from data
## the function issues a warning because some campaigns have
## not enough data points
#g.models <- gpp.bulk(flux ~ PAR + timestamp + kind, amd, amd$campaign, 
#method="Falge", min.dp=5)
#tbl8(g.models)
#
### alternative approaches to acknowledge reco when fitting GPP models
## we need only fluxes based on transparent chamber measurements (NEE)
#amg <- amd[amd$kind=="T",]
## fit gpp models and predict reco from models
#g.models.a1 <- gpp.bulk(flux ~ PAR + timestamp + t.air + t.soil2 + 
#t.soil5 + t.soil10, amg, amg$campaign, method="Falge", min.dp=5, 
#Reco.m=r.models)
#tbl8(g.models.a1)
## have a look the model fits (first 10)
#par(mfrow=c(5,6))
## select only non linear fits
#sel <- sapply(g.models.a1, function(x) class(x$mod$mg)=="nls")
#lapply(g.models.a1[sel][1:10], function(x) plot(x$mod, single.pane=FALSE))
#
## fit gpp models with providing reco data
## to do so, rerun budget.reco with other start and end points
#set.back <- data.frame(timestamp = c("2009-09-01 00:30", "2011-12-31 23:30"), 
#value = c(-999, -9999))
#set.back$timestamp <- strptime(set.back$timestamp, format="\%Y-\%m-\%d \%H:\%M")
#r.bdgt.a2 <- budget.reco(r.models, amc, set.back)
## now fit the models
#g.models.a2 <- gpp.bulk(flux ~ PAR + timestamp, amg, amg$campaign, 
#method="Falge", units = "30mins", min.dp=5, Reco.m=r.bdgt.a2$reco.flux, 
#ts.Reco = r.bdgt.a2$timestamp)
#tbl8(g.models.a2)
#
## End not run ##

}

\keyword{ univar }
