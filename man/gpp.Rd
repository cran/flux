\encoding{UTF-8}
\name{gpp}
\alias{gpp}
\title{
Model GPP from CO2 closed chamber flux data
}
\description{
Model GPP from \eqn{CO_2}{CO2} closed chamber flux data under consideration of ecosystem respiration.
}

\usage{
gpp(NEE, PAR, PAR.Temp, Reco.m, ts.NEE = NULL, ts.Reco = NULL, method = 
"Michaelis-Menten", units = "30mins", allow.offset = FALSE, virtual = FALSE, 
start.par = max(PAR), ...)
}

\arguments{
  \item{NEE}{
Numeric vector with \eqn{CO_2}{CO2} from measurements of \eqn{CO_2}{CO2} net ecosystem exchange (NEE).
}
  \item{PAR}{
Numeric vector of mean irradiation during \eqn{CO_2}{CO2} flux measurements. Can be supplied as photosynthetically active radiation (PAR) or photosynthetic photon flux density (PPFD).
}
  \item{PAR.Temp}{
Numeric vector of mean recorded temperature readings during \eqn{CO_2}{CO2} flux measurements. Has to correspond with the particular temperature used in \eqn{R_{eco}}{R_eco} modeling (e.g., air, soil in -2cm/-5cm/-10cm depth).
}
  \item{Reco.m}{
Model structure obtained from running \code{\link{reco}} or vector with estimated Reco values. The latter has to contain (at least) the values that are valid at the times of the NEE measurements. Typically the data will result from estimating hourly or half-hourly \eqn{R_{eco}}{R_eco} values using \eqn{R_{eco}}{R_eco} models derived with \code{reco} and temperature data from a climate station. 
}
  \item{ts.NEE}{
POSIXlt vector holding the timestamp of the \code{NEE} values. Has to be specified if \eqn{R_{eco}}{R_eco} values instead of an \eqn{R_{eco}}{R_eco} model are given in \code{Reco.m}. NEE and Reco values are linked to each other based on their timestamps specified in \code{ts.NEE} and \code{ts.Reco}, respectively. The two timestamps can be different and the function links NEE and Reco values based on the time interval given in \code{units}.
}
  \item{ts.Reco}{
POSIXlt vector holding the timestamp of the \eqn{R_{eco}}{R_eco} values. Has to be specified if \eqn{R_{eco}}{R_eco} values instead of an \eqn{R_{eco}}{R_eco} model are given in \code{Reco.m}. The function assumes that this is the case if \code{ts.Reco != NULL}.
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
}

\details{
The function models the relationship between \eqn{CO_2}{CO2} uptake by plants (gross primary production, GPP) and irradiation using one out of 4 methods (Falge et al. 2001). Per default the Michaelis-Menten kinetic (e.g., Schmitt et al. 2010) is used. The following models can be fitted to the data:


\eqn{GPP = \frac{GPmax * alpha * PAR}{alpha * PAR + GPmax}}{GPP = (GPmax * alpha * PAR)/(alpha * PAR + GPmax)} (Michaelis-Menten)

\eqn{GPP = \frac{alpha * PAR}{1 - \frac{PAR}{2000} + \frac{alpha*PAR}{GPmax}}}{GPP = (alpha * PAR)/(1 - (PAR/2000) + (alpha*PAR/GPmax))} (Falge)

\eqn{GPP = \frac{GPmax * alpha * PAR}{\sqrt{GPmax^2 + (alpha*par)^2}}}{GPP = (alpha * PAR * GPmax)/sqrt(GPmax^2 + (alpha*PAR)^2)} (Smith)

\eqn{GPP = GPmax * (1 - e^{\frac{alpha * PAR}{GPmax})}}{GPP = GPmax * (1 - exp((alpha*PAR) / GPmax))} (Misterlich)


with \code{PAR} the incoming light (irradiation). Note, that irradiation can be given in \code{PAR} or in PPFD although the equation states \code{PAR}. \code{GPmax} and \code{alpha} are the parameters that are fitted. \code{GPmax} refers to the maximum gross primary production at saturating or optimum light whereas alpha refers to the ecosystem quantum yield and gives the starting slope of the model.

Closed chamber measurements in the field typically capture net ecosystem exchange (\code{NEE}), which is the sum of the two opposing processes ecosystem respiration (\eqn{R_{eco}}{R_eco}) and GPP. Therefore, it is necessary to subtract modeled \eqn{R_{eco}}{R_eco} from the measured \code{NEE} to obtain GPP that can be used for the modelling against irradiance.

Real \eqn{R_{eco}}{R_eco} at the time of the \code{NEE} measurement is typically unkown because dark and light measurements cannot be taken at the same spot at the same time. Therefore, \eqn{R_{eco}}{R_eco} has to be modelled based on dark chamber or nighttime measurements (see \code{\link{reco}}). For modelling GPP from \code{NEE} chamber measurements, \code{gpp} just needs measured \code{NEE}, the associated irradiance (\code{PAR}) and temperature (\code{PAR.Temp}) values and the \eqn{R_{eco}}{R_eco} model (\code{Reco.m}). The \eqn{R_{eco}}{R_eco} model can derive from a longer period of time than the \code{NEE} data, which is often better to get more reliable models.

The Michaelis Menten fit to the \code{GPP}/\code{PAR} relationship presumes that plants (at least C3 plants) do not take up \eqn{CO_2}{CO2} when there is no irradiance. However, sometimes the \eqn{R_{eco}}{R_eco} model gives quite unrealistic \eqn{R_{eco}}{R_eco} estimates for the times of NEE measurements leading to an alleged considerable uptake of \eqn{CO_2}{CO2} under no or very low light conditions. This in turn leads to unrealistic and not well fitted GPP models. Therefore, it is possible to correct the model by not allowing an offset: \code{allow.offset = FALSE} (default). The offset is determined automatically by constructing a linear model using the data points until \code{PAR} = \code{start.par} and predicting GPP at \code{PAR} = 0. The offset is then subtracted from all GPP values and is later automatically added when doing the diagnostic plots.

The start parameters for the non-linear fit (via \code{\link[stats]{nls}}) are derived from the data itself. For alpha (initial slope of the model) the slope of the linear model of GPP against \code{PAR} constructed from  the data points until \code{PAR} = \code{start.par} is used. For \code{GPmax} the mean of the five highest GPP values is taken.

It is advisable to test various configurations regarding the \eqn{R_{eco}}{R_eco} model and testing the effect of allowing the offset. The offset is not added back to the predicted GPP data but it is returned as part of the output (see value section).
}

\value{
The function returns an object of class \code{gpp}. It is a list with the following components.

\item{mg }{ The gpp model. A \code{\link[stats]{nls}} model structure. }
\item{mr }{ The Reco model used. A \code{\link[stats]{nls}} model structure. }
\item{data }{ A three entry list with data needed for the diagnostic plot containing (1) \code{PAR.Temp} – numeric vector with the \code{PAR.Temp} values specified in the function call; (2) \code{Reco} - numeric vector of corresponding \eqn{R_{eco}}{R_eco} values estimated with the \eqn{R_{eco}}{R_eco} model (\code{Reco.m}); (3) \code{offset} – Numeric value giving the offset. }

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
gpp.m1 <- gpp(tmp$CO2.flux, tmp$PAR, tmp$t.soil2, reco.m[[1]])
# check diagnostic plot
plot(gpp.m1)

# same for spot 3
tmp <- ttf4gpp[(ttf4gpp$date=="2011-05-11") & (ttf4gpp$spot==3),]
gpp.m2 <- gpp(tmp$CO2.flux, tmp$PAR, tmp$t.soil2, reco.m[[1]])
# check diagnostic plot
plot(gpp.m2)

# same with all three spots
tmp <- ttf4gpp[(ttf4gpp$date=="2011-05-11"),]
gpp.m3 <- gpp(tmp$CO2.flux, tmp$PAR, tmp$t.soil2, reco.m[[1]])
# check diagnostic plot
plot(gpp.m3)

# do a gpp model per date with all spots each
dts <- levels(ttf4gpp$date)
# exclude measurement dates with less than 5 data points
dts <- dts[summary(ttf4gpp$date)>=5]
gpp.mdls <- lapply(dts, function(x) with(ttf4gpp[ttf4gpp$date==x,], 
gpp(CO2.flux, PAR, t.soil2, reco.m[[1]])))
# plot model diagnostic plot into one device
par(mfcol=c(3,5))
lapply(gpp.mdls, plot, single.pane=FALSE)

}

\keyword{ univar }
