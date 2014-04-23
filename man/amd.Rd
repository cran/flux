\encoding{UTF-8}
\name{amd}
\alias{amd}
\docType{data}
\title{
Closed chamber fluxes from 2009 to 2011 in the Ahlenmoor peat bog, Northeast Germany 
}
\description{
CO2 exchange rates determined with closed chamber measurements and cooresponding measurements of temperatures, photosynthetically actiove radiation, and other variables in the Ahlen-Falkenberger Moor peat bog complex from 2009 to 2011 on one specific plot (3 replicates) ....
}

\usage{data(amd)}
\format{
  A data frame with 559 observations on the following 14 variables.
  \describe{
    \item{\code{timestamp}}{POSIXlt representing the date and time for the measurements}
    \item{\code{campaign}}{Numeric. Representing IDs of the measurement campaigns, i.e. data sharing on campaign were acquired within short time period (typically one day but sometimes also two consecutive days)}
    \item{\code{plot}}{Numeric. Representing the field plot numbers}
    \item{\code{kind}}{Character vector giving the kind of chamber measurements. Either "D" for dark (opaque) chamber measurements (i.e., \eqn{R_{eco}}{R_eco} measurements) or "T" for transparent chamber measurements (i.e., NEE measurements)}
    \item{\code{flux}}{Numeric. The estimated flux rate. CO2 exchange in micromole per sqm and second}
    \item{\code{PAR}}{Numeric. Average photosynthetically active radiation during the flux measurement. Actually represents PPFD (photon flux density) in micromole per sqm and second}
    \item{\code{t.air}}{Numeric. Average air temperature during chamber measurement in °C}
    \item{\code{t.soil2}}{Numeric. Average soil temperature in 2cm depth during chamber measurement in °C}
    \item{\code{t.soil5}}{Numeric. Average soil temperature in 5cm depth during chamber measurement in °C}
    \item{\code{t.soil10}}{Numeric. Average soil temperature in 10cm depth during chamber measurement in °C}
    \item{\code{n.meas}}{Numeric giving the number of concentration measurements that were available to estimate the flux}
    \item{\code{duration}}{Factor with 282 levels giving the duration of measurement in the field, format is "\%H:\%M"}
    \item{\code{r.squared}}{Numeric. The R2s of the linear regression models that were fit to the concentration data to estimate the fluxes (with \code{\link{fluxx}})}
    \item{\code{sigmif}}{Numeric. The significance of the linear regression models}
  }
}

\source{
Saschas Doktorarbeit
}

\references{
Beetz S (2014) ???
}

\examples{
data(amd)
}

\keyword{datasets}
