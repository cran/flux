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
flux.calib(dat, columns, calib, format = "\%Y-\%m-\%d \%H:\%M:\%S", 
window = 3, buffer = 1000, n.cg = 4, rl.backup = 20, attach = FALSE)
}

\arguments{
  \item{dat}{
Object returned by \code{\link{chop}} containing gas concentration measurements for several chamber measurements.
}
  \item{columns}{
Character vector giving the names of the two columns that shall be taken from \code{dat} and from \code{calib} for extracting the calibration measurements. Typically one date and one concentration column. This also implies that they carry the same names in \code{dat} and \code{calib}. 
}
  \item{calib}{
\code{\link{data.frame}} with concentration measurements of calibration gases that have been carried out at least in part during the time the concentration measurements in \code{dat} have been achieved.
}
  \item{format}{
Character string specifying the format of dates in \code{dat} as well as in \code{calib}. Internally dates are converted to a date format R can handle (see \code{\link[base]{strptime}} for details and format options).
}
  \item{window}{
Integer value. Hours. Window around the date and time (if available) of measurement of the field greenhouse gas concentrations at the measurement device (e.g. a GC) that shall be considered for the inclusion of calibration gas measurements. If no times are given \code{window}s = 48 includes the day after the measurement date and \code{window}s >48 include the day before and the day after the measurement date.
}
  \item{buffer}{
Numeric. Concentration buffer around the range of concentration measurements in \code{dat} in which the function searches for calibration gas measurements. Defaults to 1000 (ppm or ppb, depends on gas). When only the closest calibration gas concentration shall be considered one can decrease the buffer. When real concentrations are far different from available calibration gas concentrations one might need to increase the buffer to have enough data.
}
  \item{n.cg}{
Integer. Number of calibration gas concentrations in \code{calib}.
}
  \item{rl.backup}{
Numeric value. Range limit backup value that is used in situations where no range limit can be derived from the calibration measurements. See details. Defaults to a quite reasonable 20. Deprecated.
}
  \item{attach}{
Logical. If TRUE the range limits are attached to the original data.
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
\code{\link{chop}}, \code{\link{flux}}
}

\examples{
## load example data
data(tt.pre)

## extract field concentration measurements
gcd <- tt.pre[tt.pre$sampletype_a=="P",]

## partition the data into data tables per chamber measurement
gcd.parts <- chop(gcd, factors = c("date", "spot", "veg"), 
nmes = c("date", "veg", "spot"))

## calculate range limits according to the data and the accompanying 
## calibration gas measurements
# extract and prepare calibration measurements
cal <- tt.pre[tt.pre$sampletype_a=="E",c("date_gc", "CH4ppb", "CH4Code", 
"CO2ppm", "CO2Code", "N2Oppb", "N2OCode")]
names(cal)[1] <- "date"
cal$date <- "2011-03-16"
# calculate the range limits per gas (makes no real sense with such
# a small dataset).
# CH4 range limits 
CH4.lims <- flux.calib(gcd.parts, columns = c("date", "CH4ppb"), 
calib = cal, format="\%Y-\%m-\%d", window=48, attach=FALSE, buffer=1100)
# N2O range limits 
N2O.lims <- flux.calib(gcd.parts, columns = c("date", "N2Oppb"), 
calib = cal, format="\%Y-\%m-\%d", window=48, attach=FALSE, buffer=1100)
# CO2 range limits 
CO2.lims <- flux.calib(gcd.parts, columns = c("date", "CO2ppm"), 
calib = cal, format="\%Y-\%m-\%d", window=48, attach=FALSE, buffer=1100)

## attach the range limits to the original data
gcd.parts.cal <- flux.calib(gcd.parts, columns = c("date", "CH4ppb"), 
calib = cal, format = "\%Y-\%m-\%d", attach = TRUE, window=48, buffer=1100)

}

\keyword{ manip }
\keyword{ univar }
