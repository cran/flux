\encoding{UTF-8}
\name{tt.pre}
\alias{tt.pre}
\docType{data}
\title{
One day data from dynamic closed chamber measurements in the Trebeltal
}
\description{
The data comes from the Trebeltal / Northeastern Germany and has been recorded with flexible dynamic closed chambers in March 2011. It contains concentration data from 18 chamber measurements including calibration gas measurements that have been carried out alternatingly on the GC.
}

\usage{data(tt.pre)}

\format{
  A data frame with 118 observations on the following 17 variables.
  
  \describe{
    \item{\code{year}}{numeric vector giving the year of measurement}
    \item{\code{date}}{factor giving the date of field sampling, format is "\%Y-\%m-\%d"}
    \item{\code{time}}{factor giving the time of measurement in the field, format is "\%H:\%M"}
    \item{\code{veg}}{factor with levels \code{c} \code{p} \code{t}}
    \item{\code{spot}}{numeric vector, but it is a factor giving the number of the field measurement location. The combination of \code{veg} and \code{spot} uniquely identifies the measurement locations in the site}
    \item{\code{time_min}}{numeric vector, time in minutes during the chamber placement. starts with 0 from placing the chamber}
    \item{\code{sampletype_a}}{factor with levels \code{E} \code{P} determining whether its a field concentration measurement or a calibration gas measurement}
    \item{\code{temp_dC}}{numeric vector, air temperature within chamber during measurements, taken at the same times as the concentration samples}
    \item{\code{cham_vol}}{numeric vector, chamber volume per chamber placement. Varies from chamber placement to chamber placement depending on the chamber used}
    \item{\code{cham_area}}{numeric vector giving the chamber area}
    \item{\code{date_gc}}{factor giving the date of the gc measurement, format is "\%Y-\%m-\%d"}
    \item{\code{CO2Code}}{numeric vector, quality parameter from the GC}
    \item{\code{CO2ppm}}{numeric vector, concentration of CH4 in air sample / calibration gas sample}
    \item{\code{N2OCode}}{numeric vector, quality parameter from the GC}
    \item{\code{N2Oppb}}{numeric vector, concentration of N2O in air sample / calibration gas sample}
    \item{\code{CH4Code}}{numeric vector, quality parameter from the GC}
    \item{\code{CH4ppb}}{numeric vector, concentration of CO2 in air sample / calibration gas sample}
  }
}

\details{
The 18 chamber measurements are carried out on three vegetation types (Phragmites, Typha, Carex) with two treatments each (mown / unmown). At time of measurement no treatment has been applied yet. 
}

\source{
unpublished preliminary data
}

\examples{
## load data
data(tt.pre)
## see their structure
str(tt.pre)
}

\keyword{datasets}
