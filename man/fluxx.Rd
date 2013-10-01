\encoding{UTF-8}
\name{fluxx}
\alias{fluxx}
\alias{mf.flux}
\title{
Estimate (ghg) flux rates from online dynamic closed chamber measurements in through-flow mode.
}

\description{
(Bulk) estimates of (ghg) fluxes from online concentration measurements with non-steady-state closed chambers. The function tries to find stable linear conditions in concentration change by fitting many regressions to the data and automatically detects and excludes rapid concentration fluctations.
}

\usage{
fluxx(x, var.par, subset, asterisks = FALSE, loop = "auto", ...)

mf.flux(x, var.par, method = "r2", time.unit = "S", all.through = TRUE, iv = 1, 
wndw = 0.1, pdk = 0.5, min.dp = 20, nrmse.lim = 0.1, r2.qual = 0.9, 
range.lim = 5, out.unit = "auto", elementar = FALSE, 
hardflag = list(range = TRUE), consecutive = FALSE)
}

\arguments{
  \item{x}{
A list of data tables as returned by \code{chop} or alternatively one data table (for \code{mf.flux} which will rarely be called directly by the user). Each table contains data for one chamber measurement. Minimum requirements are the columns specified in \code{var.par}.
}
  \item{var.par}{
A named list specifying the variables and parameters that are used in the estimation process and variables that should be handed through the function so that they are easily available for further analysis. Some of the names are obligatory (e.g. \code{time}, \code{volume}, \code{area}, \code{t.air}, \code{p.air}, and two or more that specify the gas column and the gas quality column), others are optional. See details and examples.
}
  \item{subset}{
An optional vector specifying a subset of concentration measurements to be used in the estimation process.
}
  \item{asterisks}{
Logical. If TRUE p-values are given as asterisks and other symbols (p<.001 = "***", .001<p<.01 = "**", .01<p<.05 = "*", .05<p<.1 = ".", p>=.1 == " ").
}
  \item{loop}{
Can be \code{TRUE}, \code{FALSE} or "auto". Determines how bulk flux estimations are done. If \code{TRUE} a \code{for}-loop is used, if \code{FALSE} \code{lapply} is used, and if "auto" the approach is switched automatically depending on the number of data tables: If \code{x} contains more than 100 data tables (chamber placements) the approach is switched from \code{lapply} to \code{for}-loop because \code{lapply} may be slow on large \code{x}.
}
  \item{\dots}{
Further arguments passed to \code{mf.flux}.
}
  \item{method}{
Character string specifying the statistic used for finding the linear part. Partial match to \code{"r2"}, \code{"rmse"}, \code{"AIC"}. Defaults to \code{"r2"}. See Details.
}
  \item{time.unit}{
Single character giving the appropriate unit of time elapsed between two concentration measurements. Will typically be seconds, thus default is \code{"S"}. Other options are \code{"M"} for minutes and \code{"H"} for hours. ATTENTION: Setting the time unit to the wrong value will result in incorrect fluxes.
}
  \item{all.through}{
Logical. When \code{TRUE}, all data columns in \code{x} other than the ones needed for flux calculation are also handed through the function so that they can be used in later steps of analysis. You may also specify one or several columns that are handed through using \code{var.par}. The \code{all.through} setting overruns the \code{var.par} settings.
}
  \item{iv}{
Numeric. Sometimes there is no time information at all but the rows in \code{x} are just numbered consecutively. The correct temporal spread is calculated inside the function when the measurement interval is specified here. Defaults to 1 which expects times to be correctly given in \code{x}.
}
  \item{wndw}{
Numeric between 0 and 1. Relative width of a moving window in which the standard deviation of the concentrations is calculated to identify high frequency fluctuations. See details and next.
}
  \item{pdk}{
Numeric between 0 and 1. Minimum proportion of data points to be kept. See details. In case one single concentration value occurs more than \code{pdk * n} times in the data (may happen under zero to very low flux conditions), all other data is assumed to represent high frequency fluctuations and flux is set to zero.
}
  \item{min.dp}{
Numeric. The minimum number of data points. Defaults to 20. If there are less rows the estimation is run anyway but a warning is issued and \code{min.dp} is automatically adjustet to n-1.
}
  \item{nrmse.lim}{
The maximum acceptable normalized root mean square error. Numeric value between 0 and 1. Defaults to 0.1. If the final best solution has a higher nrmse it is flagged accordingly.
}
  \item{r2.qual}{
Numeric between 0 and 1. Quality parameter for the model fit. The minimum acceptable \eqn{R^2}{R2} of the best fitted model. Defaults to 0.8. When the value is below quality setting a quality flag is reported.
}
  \item{range.lim}{
Numeric. The minimum range of the concentration measurements during one chamber placement. The acceptable range limit depends on the accuracy of the concentration measurements. When the range of the concentration measurements is smaller than the repeatability range of the measurement device one cannot tell real increase in concentration from random fluctuation. Therefore, if the range of the concentration measurements during one chamber placement is < \code{range.lim}, the range.flag is set to FALSE (0). See details.
}
  \item{out.unit}{
Character string determining the output unit of the flux rate mass part. The default "auto" tries to find a unit that ranges the output value between 0.01 and 10. Possible output units are "ng", "mug", "mg", or "g". "mug" stands for "\eqn{\mu}g" because non-ascii characters are not allowed in functions. Beware of varying mass units in your output when running in auto mode.
}
  \item{elementar}{
When the fluxes are wanted as element values set \code{elementar = TRUE}. Defaults to FALSE.
}
  \item{hardflag}{
Named list that controls which of the quality flags are to be hard flagged (the value is changed according to the quality flag). Only \code{range.lim} is hard flagged by default. So when the range of concentration values for a chamber measurement is smaller than the set range limit a zero flux is returned. When a flux estimation does not meet the quality requirements of any other hard-flagged quality parameter fluxx returns NA. For changing a quality flag to a hard flag just provide its name (without quotation marks) and set it to TRUE. Possible parameters to choose as hard or soft flag are \code{nrmse}, \code{range}, and \code{r2}. Further the number of measurements below ambient (\code{nomba}) can be hardflagged by setting \code{nomba = [0…nc]} with \code{nc} = number of concentration measurements for one chamber placement.
}
  \item{consecutive}{
Shall the most linear part be found by a consecutive approach starting at the first concentration reading. As soon as a stable flux is detected, it is stored. Strictly experimental.
}

}
\details{
The function is similar to \code{\link{flux}} but uses a different algorithm to identify the most linear part of the concentration development. First high frequency fluctations are omitted. Then all possible \code{pdk} * n : n consecutive concentration measurements are regressed against the corresponding times. The model with the highest r2 is chosen.

\code{var.par} specifies the variables within \code{x}  and  fixed parameters for all chamber placements that are used for the flux estimations. For obligatory \code{var.par} items see \code{\link{flux}} and examples. In contrast to \code{\link{flux}} there is just one workhorse function doing the actual estimation (\code{mf.flux}) per data table. Especially when there are many data tables in \code{x} and/or many data points per data table it takes some time. Progress is shown in the console. Each dot represents one finalized data table.
}

\value{
\code{fluxx} returns a complex object of class \code{fluxxes} that is a 2 entry list. When the object is printed to the console only the second entry is displayed in a modified form that is meant to maximize information display with small footprint for easy inspection. A table is printed to the console with three columns per gas. The first contains the quality flags (e.g. 111.9). The order is: \code{r2.f}, \code{range.f},\code{nrmse.f}, \code{nomba.f}. The first three are considered more important, and if they are '1' everything is fine. The last number digit following the full stop gives the number of concentration readings below ambient.

The \code{\link[base]{data.frame}} with the estimated flux rates contains all data needed for further analysis. The columns represent the entries in \code{fluss} of the single chamber measurements (including quality flags, see below) plus naming information according to the settings in the \code{nmes} argument of \code{\link{chop}}. \code{\link{export}} provides a simple way to export the results.

The first entry is itself a list of lists and data tables. It is called \code{flux.res}. The only one first level entry in this list contains the information for one gas which is itself a list. In this list each first level entry contains the information for one chamber measurement. It is named according to the \code{nmes}-setting in \code{\link{chop}} and contains the elements \code{fluss} (which is itself a list with the elements given below), \code{mod}, \code{out} (a list with hand through data, list items according to columns in \code{x} that have been handed trough via \code{all.through} or \code{var.par}), and \code{inn} - a data.frame with the input data that were relevant for estimating the flux (the obligatory part of \code{var.par}).

The elements of \code{fluss}:

\item{ghg }{Character. The gas species for which the flux has been estimated.}
\item{flux }{Numeric. Calculated flux rate in mass unit per m2 and hour.}
\item{r2 }{r2 of the best fitted model that has been used for flux caclulation.}
\item{nrmse }{nrmse of the best fitted model that has been used for flux calculation.}
\item{r2.f }{Logical. r2 quality flag telling whether the r2 quality setting given in \code{r2.qual} is fulfilled.}
\item{range.f }{Logical. Range quality flag telling whether the range of the concentration measurements exceeded the quality range of the measurement device that has been specified in \code{range.lim}.}
\item{nrmse.f }{Logical. nrmse quality flag telling whether the nrmse quality setting given in \code{nrmse.lim} is fulfilled (i.e. if the nrmse of the best model <= \code{nrmse.lim}).}
\item{nomba.f }{Integer. Reports the \strong{n}umber \strong{o}f \strong{m}easurements \strong{b}elow \strong{a}mbient. When one observes concentrations below ambient that might make the measurements unstable, it is possible to filter the results later and allow only a maximum acceptable number of measurements below ambient. The ambient concentration is build into the function with data from Mace Head Ireland (N2O, CH4) and global average (CO2) obtained from http://cdiac.ornl.gov/pns/current_ghg.html as of August 1st, 2011.}
\item{unit }{ The mass unit assigned. }
\item{podpu }{ Proportion (expressed as a number between 0 and 1) of data points used for constructing the linear model for estimating the flux rate. The higher the less disturbed the measurements. }

}

\references{
Nakano T (2004) A comparison of regression methods for estimating 
soil-atmosphere diffusion gas fluxes by a closed-chamber technique. 
Soil Biology and Biochemistry 36: 107-113.

Forbrich I, Kutzbach L, Hormann A, Wilmking M (2010) A comparison of 
linear and exponential regression for estimating diffusive CH4 fluxes by 
closed-chambers in peatlands. Soil Biology and Biochemistry 42: 507-515.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de}
}

\seealso{
\code{\link{gpp}} and \code{\link{reco}} for further processing of the results.
}

\examples{
\dontrun{
## load data
data(tt.nee)

## prepare flux estimation
# make parts with chop
tt.parts <- chop(tt.nee, factors=c("session", "spot"), 
nmes=c("spot", "date", "session"), min.cm=40)
# prepare var.par list (like with flux)
vp <- list(CO2 = "NEE", time = "datetime", area = "area", 
volume = "volume", t.air = "t.cham", p.air = 101325)

## do the flux estimation
# run fluxx. with lots of data it may take a while 
# (approx. 10 sec per chamber)
tt.flux <- fluxx(tt.parts, subset=c(1:30), vp, pdk=0.5, 
range.lim=3, out.unit="mg")
# inspect results table
tt.flux
# plot diagnostic plots
plot(tt.flux, dims=c(4,4), subs="spot")
# run fluxx with alternative method
tt.fluxa <- fluxx(tt.parts, subset=c(1:30), vp, pdk=0.5, 
range.lim=3, out.unit="mg", method="rmse")
# inspect results
tt.fluxa
}
}

\keyword{ univar }
