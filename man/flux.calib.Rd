\encoding{UTF-8}
\name{flux.calib}
\alias{flux.calib}
\title{
Determine calibration measurement ranges according to the dates of real measurements
}
\description{
The function basically takes calibration gas measurements and extracts the calibration gas measurements that have been carried out temporally close to a real data measurement and calculates the standard deviation of the calibration gas measurements. The obtained range limits can be used in \code{\link{flux}} as a quality parameter (via \code{range.lim}).
}

\usage{
flux.calib(dat, columns, calib, attach = TRUE, format = "\%d.\%m.\%Y", 
rl.backup = 20, calib.max = 12000, range.ext = 200, window = 7, 
calib.gas.defaults = c(300, 1000, 2000, 6000, 10000))
}

\arguments{
  \item{dat}{
Object returned by \code{\link{conz.prep}} containing gas concentration measurements for several chamber measurements.
}
  \item{columns}{
Character vector giving the names of the two columns that shall be taken from \code{dat} and from \code{calib} for extracting the calibration measurements. Typically one date and one concentration column. This also implies that they carry the same names in \code{dat} and \code{calib}. 
}
  \item{calib}{
\code{\link{data.frame}} with concentration measurements of calibration gases that have been carried out at least in part during the time the concentration measurements in \code{dat} have been achieved.
}
  \item{attach}{
Logical. If TRUE the range limits are attached to the original data.
}
  \item{format}{
Character string specifying the format of dates in \code{dat} as well as in \code{calib}. Internally dates are converted to a date format R can handle (see \code{\link[base]{strptime}} for details and format options).
}
  \item{rl.backup}{
Numeric value. Range limit backup value that is used in situations where no range limit can be derived from the calibration measurements. See details. Defaults to a quite reasonable 20.
}
  \item{calib.max}{
Numeric value. Maximum accepted concentration measurement for a calibration gas. All calibration gas measurements are filtered and any values above this limit are skipped. 
}
  \item{range.ext}{
Numeric value. Used to extent the range of the concentration measurements per chamber placement. Therefore it is in the same units as the concentration measurements themselves (typicall ppm for CO2 and ppb for CH4 and N2O).
}
  \item{window}{
Integer value. Window around the date of measurement of the field greenhouse gas concentrations at the measurement device (e.g. a GC) that shall be considered for the inclusion of calibration gas measurements.
}
  \item{calib.gas.defaults}{
Numeric or integer vector giving the pre-configured concentrations of the respective gas in the calibration gases. What really comes from the calibration gas bottles will presumably change according to fill level and other external parameters. Therefore the values can be approximations. They are just used to identify the single species of the calibration gases automatically from the achieved calibration gas concentration readings.
}
}

\details{
The function automatically detects the single species of calibration gases that have been measured. It calculates the standard deviations of the measurements per calibration gas species and than gives back an average of the calculated range limit values if there are more than one calibration gas concentrations covered by the range within the field concentration measurements per chamber placement. However, this is rather academic because a chamber measurement for which concentrations develop over the range of two or more calibration concentrations will typically not have a range limit problem.

In its actual form it is possible that there are no valid calibration measurements found for certain chamber data because the range of the chamber data (even with range extension) does not cover any of the calibration gas concentrations. In this case, the minimum range limit is assigned if \code{rl.backup} = NULL.
}

\value{
Returns a named vector with the range limits of the measurement device (as needed within \code{\link{flux}}) per chamber measurement or attaches the range limits to the original data tables that are in \code{x} and returns the altered \code{x}.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\seealso{
\code{\link{conz.prep}}, \code{\link{flux}}
}

\examples{
## load example data
data(tt.pre)

## extract field concentration measurements
gcd <- tt.pre[tt.pre$sampletype_a=="P",]

## partition the data into data tables per chamber measurement
# first tell which columns of the data you want to carry on
columns <- c("date", "veg", "spot", "time_min", "temp_dC", "cham_area", 
"cham_vol", "CO2ppm", "CO2Code", "CH4ppb", "CH4Code", "N2Oppb", "N2OCode")
# then do the partitioning
gcd.parts <- conz.prep(gcd, columns, factors = c("date", "spot", "veg"), 
nmes = c("date", "veg", "spot"))

## calculate range limits according to the data and the accompanying 
## calibration gas measurements
# extract and prepare calibration measurements
cal <- tt.pre[tt.pre$sampletype_a=="E",c("date_gc", "CH4ppb", "CH4Code", 
"CO2ppm", "CO2Code", "N2Oppb", "N2OCode")]
names(cal)[1] <- "date"
cal$date <- "2011-03-16"
# calculate the range limits per gas (makes no real sense with such
# a small dataset.
# CH4 range limits 
CH4.lims <- flux.calib(gcd.parts, columns = c("date", "CH4ppb"), 
calib = cal, format="\%Y-\%m-\%d")
# N2O range limits 
N2O.lims <- flux.calib(gcd.parts, columns = c("date", "N2Oppb"), 
calib = cal, format="\%Y-\%m-\%d")
# CO2 range limits 
CO2.lims <- flux.calib(gcd.parts, columns = c("date", "CO2ppm"), 
calib = cal, format="\%Y-\%m-\%d")

## attach the range limits to the original data
gcd.parts.cal <- flux.calib(gcd.parts, columns = c("date", "CH4ppb"), 
calib = cal, format = "\%Y-\%m-\%d", attach = TRUE)

}

\keyword{ manip }
\keyword{ univar }
