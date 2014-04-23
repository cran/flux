\encoding{UTF-8}
\name{amc}
\alias{amc}
\docType{data}
\title{
Climate station data from 2009 to 2011 in the Ahlenmoor peat bog, Northeast Germany 
}

\description{
Climatic variables measured from 2009 to 2011 in the Ahlenmoor peat bog, Northeast Germany as part of a closed chamber measurement study on GHG exchange.
}

\usage{data(amc)}
\format{
  A data frame with 43197 observations on the following 8 variables.
  \describe{
    \item{\code{date}}{Factor giving the date of field sampling, format is "\%Y-\%m-\%d".}
    \item{\code{time}}{Factor giving the time of measurement in the field, format is "\%H:\%M:\%S".}
    \item{\code{t.air}}{Numeric. Average air temperatures in °C}
    \item{\code{t.soil2}}{Numeric. Average soil temperatures in 2cm depth in °C}
    \item{\code{t.soil5}}{Numeric. Average soil temperatures in 5cm depth in °C}
    \item{\code{t.soil10}}{Numeric. Average soil temperatures in 10cm depth in °C}
    \item{\code{PAR}}{Numeric. Average half hourly photosynthetically active radiation during the flux measurement. Actually represents PPFD (photon flux density) in micromole per sqm and second}
    \item{\code{timestamp}}{POSIXlt representing the date and time for the measurements}
  }
}

\source{
Diss Sascha
}

\references{
Beetz (2014) ...
}

\examples{
data(amc)
}

\keyword{datasets}
