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
gpp(NEE, PAR, PAR.Temp, Reco.m, allow.offset = FALSE)
}

\arguments{
  \item{NEE}{
Numeric vector with \eqn{CO_2}{CO2} fluxes obtained under light conditions.
}
  \item{PAR}{
Numeric vector with PAR readings that where obtained during \eqn{CO_2}{CO2} flux measurements.
}
  \item{PAR.Temp}{
Numeric vector with temperature readings that where obtained during \eqn{CO_2}{CO2} flux measurements. The relevant temperature (e.g., air, soil in -2cm/-5cm/-10cm depth) depends on the ecosystem under study.
}
  \item{Reco.m}{
Model structure obtained from running \code{\link{reco}}
}
  \item{allow.offset}{
Logical. Shall GPP values other than 0 be allowed at \code{PAR} = 0? See details.
}
}

\details{
There is a well known relationship between \eqn{CO_2}{CO2} uptake by plants and photosynthetically active radiation that can be modelled using the Michaelis-Menten kinetic (e.g., Schmitt et al. 2010). This is what the function does. However, in the field net ecosystem exchange (NEE) is measured. NEE largely consists of two components: Ecosystem respiration (\eqn{R_{eco}}{R_eco}) and the uptake of \eqn{CO_2}{CO2} by plants (gross primary production, GPP). Therefore, it is necessary to subtract \eqn{R_{eco}}{R_eco} from the measured NEE to obtain GPP that can be used for the modelling against PAR.

Real \eqn{R_{eco}}{R_eco} at the time of the NEE measurement is typically unkown because neither with closed chambers nor with eddy covariance is it possible to do dark and light measurements at the same time. Therefore, \eqn{R_{eco}}{R_eco} has to be modelled based on dark chamber or nighttime measurements (see (\code{\link{reco}}). For modelling GPP from NEE chamber measurements with \code{gpp}, the function just needs the measured \code{NEE}, the associated \code{PAR} and \code{PAR.Temp}erature values and the \eqn{R_{eco}}{R_eco} model (\code{Reco.m}) that shall be used for estimating \eqn{R_{eco}}{R_eco} values using the \code{PAR.Temp} data. It doesn't matter whether the \eqn{R_{eco}}{R_eco} model covers the same period as the \code{NEE} data or is derived from a longer period of time. The latter is often better to get more reliable models.

The Michaelis Menten fit to the \code{GPP}/\code{PAR} relationship bases on the idea that plants (at least C3 plants) do not take up \eqn{CO_2}{CO2} when there is no \code{PAR}. However, sometimes it happens that the Reco model gives quite unrealistic estimates for \eqn{R_{eco}}{R_eco} during NEE measurements leading to an alleged considerable uptake of \eqn{CO_2}{CO2} under no or very low light conditions. This in turn leeds to unrealistic and not well fitted models. Therefore, it is possible to correct the offset (default) by not allowing an offset: \code{allow.offset = FALSE}. The offset is determined automatically by constructing a linear model using the data points until \code{PAR} = 500 and predicting GPP at \code{PAR} = 0. The offset is then subtracted from all GPP values and is later automatically added when doing the diagnostic plots. 

The start parameters for the non-linear fit (via \code{\link[stats]{nls}}) are derived from the data itself. For alpha (initial slope of the curve) the slope of  a linear model of GPP against \code{PAR} constructed from  the data points until \code{PAR} = 500 is used. For GPmax the mean of the five highest GPP values is taken.

It is advisable to test various configurations regarding the \eqn{R_{eco}}{R_eco} model and testing the effect of allowing the offset.  
}

\value{
The function returns an object of class \code{gpp}. It is a list with the following components.

\item{mg }{ The gpp model. A \code{\link[stats]{nls}} model structure. }
\item{mr }{ The Reco model used. A \code{\link[stats]{nls}} model structure. }
\item{data }{ A three entry list with data needed for the diagnostic plot containing (1) \code{PAR.Temp} – numeric vector with the \code{PAR.Temp} values specified in the function call; (2) \code{Reco} - numeric vector of corresponding \eqn{R_{eco}}{R_eco} values estimated with the \eqn{R_{eco}}{R_eco} model (\code{Reco.m}); (3) \code{offset} – Numeric value giving the offset. }

}

\references{
Schmitt M, Bahn M, Wohlfahrt G, Tappeiner U, Cernusca A. 2010. Land use affects the net ecosystem CO2 exchange and its components in mountain grasslands. Biogeosciences, 7:2297-2309.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de}
}

\seealso{
\code{\link{reco}}
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
# exclude measurement dates with less then 5 data points
dts <- dts[summary(ttf4gpp$date)>=5]
gpp.mdls <- lapply(dts, function(x) with(ttf4gpp[ttf4gpp$date==x,], 
gpp(CO2.flux, PAR, t.soil2, reco.m[[1]])))
# plot model diagnostic plot into one device
par(mfcol=c(3,5))
lapply(gpp.mdls, plot, single.pane=FALSE)

}

\keyword{ univar }
