\encoding{UTF-8}
\name{flux}
\alias{flux}
\alias{flux.odae}
\alias{flux.conv}
\title{
Estimate gas flux rates using dynamic closed chamber data
}

\description{
The functions are used to prepare and carry out gas flux rate estimation with dynamic closed chamber data for one gas at a time. Either you can do the calculation in two steps or in just one single step. \code{flux} is a convenience wrapper for the other two functions that should be the best for most users. 
}

\usage{
flux(x, columns, M = "CH4", min.allowed = 3, max.nrmse = 0.1, 
nrmse.lim = 0.2, range.lim = 30, rl.backup = 30, r2.qual = 0.8, 
nomba = NULL, in.unit = "ppb", out.unit = "mg")

flux.odae(dat, min.allowed = 3, max.nrmse = 0.1, start = 1)

flux.conv(fl.dat, ch.area, M = "CH4", range.lim = 30, r2.qual = 0.8, 
nrmse.lim = 0.2, nomba = NULL, in.unit = "ppm", out.unit = "mg")
}

\arguments{
  \item{x}{
A list with data tables as returned by \code{conz.prep}. Each table contains data for one chamber mesurement. Usually there will be identifying columns, date, time and time sequence columns, and columns containing the concentrations of the greenhouse gases for each time step. At the minimum the columns specified in the argument \code{columns} have to be there.
}
  \item{columns}{
The columns that are used in the calculation process. They should represent in this sequence: gas concentration, closing time (of the chamber) in minutes, gas concentration quality flag (might be provided by your measurement device, if not just make it 0 but provide the column anyway), chamber volume (constant throughout one chamber placement), temperature within the chamber headspace during measurements (may change during chamber placement, an average is calculated within the function), chamber area (constant during chamber placement). See example data.  
}
  \item{M}{
Character string identifying the greenhouse gas for which concentration measurements are handled. Can be \code{"CH4"}, \code{"N2O"} or \code{"CO2"}. Defaults to \code{"CH4"}. This argument determines the molar weight that is used in the flux calculation with \code{\link{gflux}} and the ambient level of the gas that is added to the diagnostic plots plotted with \code{plot.fluss} and is used for the determination of the number of measurements below ambient (\code{nomba}).
}
  \item{min.allowed}{
Minimum number of concentration measurements that are allowed for the determination of the best model fit during the outlier detection and elimination (odae) procedure. Integer between 3 and the number of concentration measurements during one chamber placement.
}
  \item{max.nrmse}{
The maximum acceptable normalized root mean square error for configurations with higher numbers of concentration measurements. Numeric value between 0 and 1. Defaults to 0.1. Above that value lower numbers of concentration measurements down to \code{min.allowed} are considered. See details.
}
  \item{nrmse.lim}{
The main quality parameter for the model fit. The maximum acceptable normalized root mean square error. Numeric between 0 and 1. Default is 0.2. If the best fit for one chamber measurement is above this value, the function still returns a flux value for that chamber measurement but reports \code{FALSE} in the \code{nrmse.check} quality flag. See details and value.
}
  \item{range.lim}{
The minimum range of the concentration measurements during one chamber placement. You can either give one numerical value or a numeric vector with the same length as \code{x} giving different range limits for each chamber placement (for instance obtained by \code{\link{flux.calib}}). The acceptable range limit depends on the accuracy of the concentration measurements. When the range of the concentration measurements is smaller than the repeatability range of the measurement device (e.g., a gas chromatograph) one cannot tell real increase in concentration from random fluctuation. Therefore, if the range of the concentration measurements during one chamber placement is < \code{range.lim}, the flux rate is declared to be 0. See details.
}
  \item{rl.backup}{
When a vector of range limits per chamber placement is obtained via \code{range.lim} it might occur that there are NA values (when there are no calibration measurements to support a regular measurement - for details see \code{\link{flux.calib}}). Therefore, you can specify a backup range limit here. Expects a numeric. Defaults to 30.
}
  \item{r2.qual}{
Alternative quality parameter for the model fit. The minimum acceptable r2 of the best fitted model. Numeric between 0 and 1. Defaults to 0.8. A chamber measurement is not defined invalid when the value is below the quality setting. Instead, a quality flag is reported (\code{TRUE} - r2-quality setting fulfilled, \code{FALSE} - r2-quality setting not fulfilled).
}
  \item{nomba}{
Quality setting: Sets the acceptable \strong{n}umber \strong{o}f \strong{m}easurements \strong{b}elow \strong{a}mbient. Defaults to the number of concentration measurements during one chamber placement. When one observes concentrations below ambient that might make the measurements unstable, it is possible to set a maximum acceptable number of measurements below ambient here. The ambient concentration is determined by the argument \code{M}. The flux rate is calculated anyway but a nomba quality flag is reported similarly to the r2-quality flag).
}
  \item{in.unit}{
Input unit of the concentration measurements. Character string. Can be "ppm" or "ppb". Defaults to "ppm".)
}
  \item{out.unit}{
Output unit of the flux rate mass part. Character string. Can be "ng", "mug", "mg", or "g". Defaults to "mg". "mug" stands for "\eqn{\mu}g" because non-ascii characters are not allowed in functions.
}
  \item{dat}{
One data table for one chamber placement. See \code{x} and \code{columns} for details.
}
  \item{start}{
Sometimes several concentration samples are taken at time 0. If so, you can specifiy here which of those shall be taken for the flux estimation.
}
  \item{fl.dat}{
An object with the same structure as returned by \code{\link{flux.odae}}. See details and value.
}
  \item{ch.area}{
Chamber area. For historical reasons this can be specified separately in the function \code{flux.conv}. When flux rate estimation is carried out using \code{flux} it is derived from the data tables.
}
}

\details{
Typically it will be most convenient to use \code{flux} on objects returned by \code{\link{conz.prep}} (i.e. on lists of data tables that contain all necessary data per chamber measurement including supporting information). \code{flux} simply wraps \code{flux.odae} and \code{flux.conv} applied on lists of chamber measurement data tables into one function. Thus, the data of a one day field campaign or a year of chamber measurements can easily handled by simply running two functions (\code{conz.prep} and \code{flux}) consecutively per gas species. See example.

If the flux estimation is carried out in two steps it will typically be carried out on a list structure as returned by \code{\link{conz.prep}}. Therefore, it is used within a \code{\link{lapply}} call. For details see examples. However, the functions \code{flux.odae} and \code{flux.conv} are  designed to be carried out on single data tables (\code{\link[base]{data.frame}}) per chamber measurement. 

First \code{flux.odae} is run. It simply tries to find the best model fit (for th etime being using linear regression) for the series of concentration measurements that are given in \code{dat}. This \code{\link[base]{data.frame}} has to consist of five columns that give (in that order): gas concentration, closing time of the chamber in minutes, gas concentration quality flag, chamber volume, temperature within the chamber headspace during measurements (may change during chamber placement). See example data. 

At the moment the optimization bases on linear regression. All possible models with n (= total number of concentration measurements per chamber placement) to \code{min.allowed} number of concentration measurements are fitted and the best fit is evaluated in a stepwise procedure. The normalized root mean square error is used as the quality criterion for the outlier detection and elimination procedure. All model fits with a nrmse <= \code{max.nrmse} are extracted and ranked according to the number of concentration measurements (decreasing) and to the nrmse (increasing). The first ranked model is stored along with the original data table and some other information. Therefore a model with e.g. a nrmse of 0.081 constructed from 5 concentration measurements wins against a model with a nrmse of 0.07 with only 4 concentration measurements. This reflects the idea that models with nrmse <= \code{max.nrmse} already represent a sufficient fit and do not have "outliers" that must be eliminated. 

In case no model has a nrmse <= \code{max.nrmse}, the models are simply ranked according to their nrmse and the model with the lowest nrmse wins and is stored. In that way outliers are detected and exluded. \code{flux.odae} returns a complex object that contains most of the necessary information for the \code{flux.conv} function and also carries information that is later needed for the plot functions (\code{plot.flux} and \code{plot.fluss}). 

The flux calculation is then carried out with the function \code{flux.conv}. It takes the object returned by \code{flux.odae} and additional information (chamber area, gas species, several quality settings and in- as well as output units) and calculates the flux rates. Further several quality checks (r2 check, range check, nrmse check, nomba check; for details see Value) are carried out and quality flags are reported along with the fluxes in the output. It is best when all quality flags are returned \code{TRUE}. Depending on the application quality requirements might vary. Therefore, we function rather report quality flags than excluding data. 
}

\value{
\code{flux} returns an object of class \code{fluss} that is a two entry list. 

The second entry of that list is a \code{\link[base]{data.frame}} with the estimated flux rates. The columns represent the entries in \code{fluss} of the single chamber measurements (see below) plus naming information according to the settings in the \code{nmes} argument of \code{\link{conz.prep}}.

The first entry is itself a list of lists and data tables. It is called \code{flux.res} and is also the object that is returned by \code{flux.conv}. Each first level entry in this list contains the information for one chamber measurement. It is named according to the \code{nmes}-setting in \code{\link{conz.prep}} and contains the elements \code{fluss} (which is itself a list with the elements given below), \code{fl.dat} (equals the object returned by \code{flux.odae}; see below), and \code{output.unit} which provides information on the output mass unit of the flux rate that is handed over to the function \code{\link{plot.fluss}}.

The elements of \code{fluss}:

\item{ghg }{Character. The gas species for which the flux has been estimated (\code{M} is just forwarded).}
\item{flux }{Numeric. Calculated flux rate in mass unit per m2 and hour.}
\item{r2.check }{Logical. r2 quality flag telling whether the r2 quality setting given in \code{r2.qual} is fulfilled.}
\item{range.check }{Logical. Range quality flag telling whether the range of the concentration measurements exceeded the quality range of the measurement device that has been specified in \code{range.lim}. In case \code{range.check = FALSE}, \code{flux} is automatically returned as 0.}
\item{nrmse.check }{Logical. nrmse quality flag telling whether the nrmse quality setting given in \code{nrmse.lim} is fulfilled (i.e. if the nrmse of the best model <= \code{nrmse.lim}).}
\item{nomba.check }{Logical. nomba quality flag telling whether the number of concentration measurements has been <= the setting in \code{nomba}.}
\item{r2 }{The actual r2 of the best fitted model that has been used for flux caclulation.}
\item{nrmse }{The actual nrmse of the best fitted model that has been used for flux caclulation.}

The elements of \code{fl.dat} that is also the object returned by \code{flux.odae} are:

\item{lm4flux }{Complex object. The best fitting model as reported by \code{\link{lm}}. It builds the basis for the calculation of the flux rate via \code{flux.conv}.}
\item{row.select }{Integer vector giving the indices of the rows of the data table that have been used to construct the best fitting model. This information is later used in the plotting functions \code{\link{plot.flux}} and \code{\link{plot.fluss}}. }
\item{orig.dat }{\code{\link[base]{data.frame}} with the original data provided according to arguments \code{x} and \code{columns}. }

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
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>, Franziska Koebsch <franziska.koebsch@uni-rostock.de>
}

\seealso{
\code{\link{conz.prep}}, \code{\link{flux.calib}}, \code{\link{gflux}}, \code{\link{plot.fluss}}
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

## calculate flux rates for methane
# first define CH4 range limit (alternatively use flux.calib)
CH4.lim <- 30
# do the flux rate estimation
flux.CH4 <- flux(gcd.parts, columns=c("CH4ppb", "time_min", "CH4Code", 
"cham_vol", "temp_dC", "cham_area"), M = "CH4")

## look at the results table
flux.CH4$flux.table

## plot the concentration-change-with-time-plots as kind of diagnostic
plot(flux.CH4, subs = NULL, dims = c(3,6))

}

\keyword{ univar }
