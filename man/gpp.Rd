\encoding{UTF-8}
\name{gpp}
\alias{gpp}
\alias{gpp2}
\title{
Model GPP from CO2 closed chamber flux data
}
\description{
Model GPP from \eqn{CO_2}{CO2} closed chamber flux data under consideration of ecosystem respiration. Four different methods are available: Providing one global Reco model, providing several Reco models, providing estimated Reco fluxes via function \code{gpp} or extracting Reco fluxes from real measurements via \code{gpp2}. Timestamps are used to assign Reco data to the respective NEE data. In the latter case they have to be provided alongside the Reco fluxes.
}

\usage{
gpp(NEE, PAR, ts.NEE, PAR.Temp, Reco.m, ts.Reco = NULL, 
method = "Michaelis-Menten", units = "30mins", allow.offset = FALSE, 
virtual = FALSE, start.par = max(PAR), ...)

gpp2(NEE, PAR, ts.NEE, oot, oot.id = c("D", "T"), 
method = "Michaelis-Menten", allow.offset = FALSE, 
virtual = FALSE, start.par = max(PAR), ...)
}

\arguments{
  \item{NEE}{
Numeric vector with \eqn{CO_2}{CO2} fluxes from measurements of \eqn{CO_2}{CO2} net ecosystem exchange (NEE).
}
  \item{PAR}{
Numeric vector of mean irradiation during \eqn{CO_2}{CO2} flux measurements. Can be supplied as photosynthetically active radiation (PAR) or photosynthetic photon flux density (PPFD).
}
  \item{ts.NEE}{
POSIXlt vector holding the timestamp of the \code{NEE} values. NEE and Reco values are linked to each other based on their timestamps specified in \code{ts.NEE} and \code{ts.Reco} or elsewhere (depending on method). The two timestamps don't have to match exactly. In that case, the function links NEE and Reco values based on the time interval given in \code{units}.
}
  \item{PAR.Temp}{
Either numeric vector of mean recorded temperature readings during \eqn{CO_2}{CO2} flux measurements or data.frame with several temperature records (if Reco.m is provided as an object resulting from running \code{\link{reco.bulk}}). In case of the latter, appropriate temperatures are extracted based on the \code{which.Temp} parameter that is stored to the model structure that is returned by \code{\link{reco.bulk}}. Therefore names have to correspond with the particular temperature variable names used in \eqn{R_{eco}}{R_eco} modeling (e.g., air, soil in -2cm/-5cm/-10cm depth). See details.
}
  \item{Reco.m}{
Model structure obtained from running \code{\link{reco.bulk}} or \code{\link{reco}} or vector with estimated Reco values. The latter has to contain (at least) the values that are valid at the times of the NEE measurements. Typically the data will result from estimating hourly or half-hourly \eqn{R_{eco}}{R_eco} values using \code{\link{budget.reco}} on objects of class \code{"breco"} derived from running \code{\link{reco.bulk}}. See details.
}
  \item{ts.Reco}{
POSIXlt vector holding the timestamp of the \eqn{R_{eco}}{R_eco} values. Has to be specified if \eqn{R_{eco}}{R_eco} values instead of an \eqn{R_{eco}}{R_eco} model are given in \code{Reco.m}. The function assumes that this is the case if \code{ts.Reco != NULL}. See details.
}
  \item{method}{
The function knows several equations to model the relationship between gpp and irradiation. At the moment \code{"Michaelis-Menten"}, \code{"Falge"}, \code{"Smith"}, and \code{"Misterlich"} are implemented which are all discussed in Falge et al. 2001. Partial matching is applied. Defaults to \code{"Michaelis-Menten"}. See details for equations.
}
  \item{units}{
Character string specifying how \code{ts.NEE} shall be rounded. If \code{Reco.m} holds values instead of an \eqn{R_{eco}}{R_eco} model structure the NEE and \eqn{R_{eco}}{R_eco} values are matched based on their timestamp after \code{ts.NEE} has been rounded according to \code{units}.
}
  \item{allow.offset}{
Logical. Shall GPP values other than 0 be allowed at zero irradiation? See details.
}
  \item{virtual}{
Logical. If \code{TRUE}, virtual NEE data are generated that show a typical saturation curve with saturation at mean \code{NEE}. Can be used in bulk gpp modeling to allow falling back to a mean model.
}
  \item{start.par}{
Numeric between 0 and \code{max(PAR)}. All data points with \code{PAR <= start.par} are used to obtain a start value for alpha via linear regression. Defaults to \code{max(PAR)}. See details.
}
  \item{\dots}{
Any arguments passed to \code{\link[stats]{nls}} which is used internally to do the model fitting.
}
  \item{oot}{
Vector of length = \code{length(NEE)} specifying which of the measured fluxes derive from opaque (\eqn{R_{eco}}{R_eco}) and which derive from transparent (NEE) chamber measurements. \code{gpp2} uses this to extract corresponding \eqn{R_{eco}}{R_eco} values for calculating GPP from NEE before fitting the models (Approach 1, see details.)
}
  \item{oot.id}{
Vector of length 2 that specifies which of the flux values derive from opaque (first value, i.e. \eqn{R_{eco}}{R_eco} measurements) and which derive from transparent (second value, i.e. NEE measurements) chamber measurements when \code{data} contains both. May be character, factor, or numeric. See details.
}
}

\details{
The function models the relationship between \eqn{CO_2}{CO2} uptake by plants (gross primary production, GPP) and irradiation using one out of 4 methods (Falge et al. 2001). Per default the Michaelis-Menten kinetic (e.g., Schmitt et al. 2010) is used. The following models can be fitted to the data:


\eqn{GPP = \frac{GPmax * alpha * PAR}{alpha * PAR + GPmax}}{GPP = (GPmax * alpha * PAR)/(alpha * PAR + GPmax)} (Michaelis-Menten)

\eqn{GPP = \frac{alpha * PAR}{1 - \frac{PAR}{2000} + \frac{alpha*PAR}{GPmax}}}{GPP = (alpha * PAR)/(1 - (PAR/2000) + (alpha*PAR/GPmax))} (Falge)

\eqn{GPP = \frac{GPmax * alpha * PAR}{\sqrt{GPmax^2 + (alpha*par)^2}}}{GPP = (alpha * PAR * GPmax)/sqrt(GPmax^2 + (alpha*PAR)^2)} (Smith)

\eqn{GPP = GPmax * (1 - e^{\frac{alpha * PAR}{GPmax})}}{GPP = GPmax * (1 - exp((alpha*PAR) / GPmax))} (Misterlich)


with \code{PAR} the incoming light (irradiation). Note, that irradiation can be given in \code{PAR} or in PPFD although the equation states \code{PAR}. \code{GPmax} and \code{alpha} are the parameters that are fitted. \code{GPmax} refers to the maximum gross primary production at saturating or optimum light whereas alpha refers to the ecosystem quantum yield and gives the starting slope of the model.

Transparent closed chamber measurements in the field typically capture net ecosystem exchange (\code{NEE}), which is the sum of the two opposing processes ecosystem respiration (\eqn{R_{eco}}{R_eco}) and GPP. Therefore, it is necessary to subtract modeled \eqn{R_{eco}}{R_eco} from the measured \code{NEE} to obtain GPP that can be used for the modelling against irradiance.

Real \eqn{R_{eco}}{R_eco} at the time of the \code{NEE} measurement is typically unkown because dark and light measurements cannot be taken at the same spot at the same time. Therefore, \eqn{R_{eco}}{R_eco} has to be modelled based on dark chamber or nighttime measurements (see \code{\link{reco}}). For modelling GPP from \code{NEE} chamber measurements, \code{gpp} just needs measured \code{NEE}, the associated irradiance (\code{PAR}) and temperature (\code{PAR.Temp}) values and the \eqn{R_{eco}}{R_eco} model(s) (\code{Reco.m}). The \eqn{R_{eco}}{R_eco} model(s) can derive from a longer period of time than the \code{NEE} data, which is often better to get more reliable models. In contrast, \code{gpp2} extracts \eqn{R_{eco}}{R_eco} fluxes from actual measurements.

Approaches to assigning \eqn{R_{eco}}{R_eco} values: 

\bold{Approach 1:} Extract corresponding \eqn{R_{eco}}{R_eco} fluxes from the provided data that are assigned to corresponding NEE values via their timestamp: For this approach \code{NEE} has to contain both NEE and \eqn{R_{eco}}{R_eco} fluxes. \code{oot} has to be specified as a vector that indicates whether the respective fluxes were measured as \code{NEE} (transparent chamber) or \code{Reco} (opaque chamber or low PAR). In addition \code{oot.id} may have to be changed accordingly. \code{\link{gpp2}} is used for fitting the models.

\bold{Approach 2:} If \code{Reco.m} is specified as a vector containing modelled \eqn{R_{eco}}{R_eco} values these are used to calculate GPP = NEE + Reco. The correct \eqn{R_{eco}}{R_eco} values are assigned to the appropriate \code{NEE} values by rounding the timestamp of the latter (given in \code{ts.NEE}) according to the time lapse of the \eqn{R_{eco}}{R_eco} values and then merging both on the respective timestamps. Therefore \code{ts.Reco} has to be specified while \code{PAR.Temp} is ignored.

\bold{Approach 3:} If just one \eqn{R_{eco}}{R_eco} model is provided as an object of class \code{"reco"} resulting from running \code{\link{reco}} this is used to predict \eqn{R_{eco}}{R_eco} at the times of the \code{NEE} measurements with the temperatures provided in \code{PAR.Temp} as \code{new.data}. \code{PAR.Temp} has to be specified as a vector of length = length(NEE). \code{ts.Reco} must not be specified.

\bold{Approach 4:} If several \eqn{R_{eco}}{R_eco} models are provided as an object of class \code{"breco"} resulting from running \code{\link{reco.bulk}} these are used to predict \eqn{R_{eco}}{R_eco} at the times of the \code{NEE} measurements with the temperatures provided in \code{PAR.Temp} as \code{new.data}. \code{PAR.Temp} has to be provided as a data.frame with all temperature variables that were used when obtaining the \eqn{R_{eco}}{R_eco} models via \code{\link{reco.bulk}} with \code{ncol(PAR.Temp)} = \code{length(NEE)}. The appropriate temperatures are assigned using the parameter \code{which.Temp} that is reported with each model in an object of class \code{"breco"}. \code{ts.Reco} must not be specified.

The Michaelis Menten fit to the \code{GPP}/\code{PAR} relationship presumes that plants (at least C3 plants) do not take up \eqn{CO_2}{CO2} when there is no irradiance. However, sometimes the \eqn{R_{eco}}{R_eco} model gives quite unrealistic \eqn{R_{eco}}{R_eco} estimates for the times of NEE measurements leading to an alleged considerable uptake of \eqn{CO_2}{CO2} under no or very low light conditions. This in turn leads to unrealistic and not well fitted GPP models. Therefore, it is possible to correct the model by not allowing an offset: \code{allow.offset = FALSE} (default). The offset is determined automatically by constructing a linear model using the data points until \code{PAR} = \code{start.par} and predicting GPP at \code{PAR} = 0. The offset is then subtracted from all GPP values and is later automatically added when doing the diagnostic plots.

The start parameters for the non-linear fit (via \code{\link[stats]{nls}}) are derived from the data itself. For alpha (initial slope of the model) the slope of the linear model of GPP against \code{PAR} constructed from  the data points until \code{PAR} = \code{start.par} is used. For \code{GPmax} the mean of the five highest GPP values is taken.

It is advisable to test various configurations regarding the \eqn{R_{eco}}{R_eco} model and testing the effect of allowing the offset. ATTENTION: The offset is not added back to the predicted GPP data but it is returned as part of the output (see value section). Therefore, if the model parameters and model formula are used to predict GPP fluxes, the offset has to be added manually.
}

\value{
The function returns an object of class \code{gpp} (for \code{ts.Reco} != \code{NULL}) or of class \code{gpp2} (for \code{ts.Reco} = \code{NULL}). It is a list with the following components.

\item{mg }{ The gpp model. A \code{\link[stats]{nls}} model structure. }
\item{mr }{ The Reco model used. A \code{\link[stats]{nls}} model structure. }
\item{data }{ Either a three entry list (with \code{ts.Reco != NULL}) or a 4 entry list (with \code{ts.Reco = NULL}) }
\item{dat}{ \code{data.frame} (see below for more). }
\item{offset }{ Numeric value giving the offset. }
\item{start }{ List with the start values for the \code{gpp} modelling. }
\item{PAR.Temp }{ Numeric vector with the \code{PAR.Temp} values specified in the function call. Only reported if \code{ts.Reco = NULL}. }

The \code{data.frame} in \code{dat} contains the following columns:

\item{NEE }{ \code{NEE} values. }
\item{GPP }{ Corresponding \code{GPP} values. }
\item{Reco }{ Corresponding \eqn{R_{eco}}{R_eco} values. }
\item{PAR }{ Corresponding \code{PAR} values. }
\item{timestamp }{ Corresponding timestamps. }
\item{mins }{ Temporal distance to next reco value. Always 0 but reported for consistency with \code{\link{gpp2}}. }
\item{Reco }{ Numeric vector of corresponding \eqn{R_{eco}}{R_eco} values estimated with the \eqn{R_{eco}}{R_eco} model (\code{Reco.m}). }

}

\references{
Falge E, Baldocchi D, Olson R, Anthoni R, et al. 2001. Gap filling strategies for defensible annual sums of net ecosystem exchange. Agricultural and Forest Meteorology, 107:43-69.

Schmitt M, Bahn M, Wohlfahrt G, Tappeiner U, Cernusca A. 2010. Land use affects the net ecosystem CO2 exchange and its components in mountain grasslands. Biogeosciences, 7:2297-2309.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de}
}

\seealso{
\code{\link{reco}}, \code{\link{fluxx}}
}

\examples{
## load data
data(tt.flux)

## make timestamp
tt.flux$timestamp <- strptime(paste(tt.flux$date, tt.flux$time), 
format="\%Y-\%m-\%d \%H:\%M:\%S")

## model reco with Arrhenius type model
# extract data and omit estimated fluxes with both the nrmse 
# and the r2 flag set to 0
ttf <- tt.flux[!(tt.flux$CO2.r2.f + tt.flux$CO2.nrmse.f) == 0, ]

# extract table with flux data for reco modeling
ttf4reco <- subset(ttf, kind > 4)

# omit CO2 fluxes below zero
ttf4reco <- ttf4reco[ttf4reco$CO2.flux >= 0,]

# plot reco data
plot(CO2.flux ~ t.air, data=ttf4reco)

# check for the best temperature for reco modelling
temps <- c("t.air", "t.soil2", "t.soil5", "t.soil10")
sapply(temps, function(x) lapply(reco(ttf4reco$CO2.flux, 
ttf4reco[,x], method="arr"), AIC))

# take the temperature in soil 2 cm
reco.m <- reco(ttf4reco$CO2.flux, ttf4reco$t.soil2, method="arr")

# inspect
reco.m

## model gpp
# extract table with flux data for gpp modeling
ttf4gpp <- subset(ttf, kind < 4)

# do a single gpp model for a measurement day using data of spot 2
tmp <- ttf4gpp[(ttf4gpp$date=="2011-05-11") & (ttf4gpp$spot==2),]
gpp.m1 <- gpp(tmp$CO2.flux, tmp$PAR, tmp$timestamp, tmp$t.soil2, 
reco.m[[1]])
# check diagnostic plot
plot(gpp.m1)

# same for spot 3
tmp <- ttf4gpp[(ttf4gpp$date=="2011-05-11") & (ttf4gpp$spot==3),]
gpp.m2 <- gpp(tmp$CO2.flux, tmp$PAR, tmp$timestamp, tmp$t.soil2, 
reco.m[[1]])
# check diagnostic plot
plot(gpp.m2)

# same with all three spots
tmp <- ttf4gpp[(ttf4gpp$date=="2011-05-11"),]
gpp.m3 <- gpp(tmp$CO2.flux, tmp$PAR, tmp$timestamp, tmp$t.soil2, 
reco.m[[1]])
# check diagnostic plot
plot(gpp.m3)

}

\keyword{ univar }
