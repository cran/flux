\encoding{UTF-8}
\name{fluxx}
\alias{fluxx}
\alias{mf.flux}
\title{
Estimate (ghg) flux rates from online dynamic closed chamber measurements in through-flow mode.
}

\description{
(Bulk) estimates of (ghg) flux rates from medium frequency (mf) online concentration measurements with dynamic closed chambers in through-flow mode. The function tries to find stable linear conditions in concentration change and automatically detects and excludes concentration fluctations.
}

\usage{
fluxx(x, var.par, subset, asterisks = FALSE, loop = "auto", ...)

mf.flux(x, var.par, time.unit = "S", all.through = TRUE, iv = 1, 
wndw = 0.1, pdk = 0.5, min.dp = 20, nrmse.lim = 0.1, r2.qual = 0.9, 
range.lim = 5, out.unit = "auto", elementar = FALSE, 
hardflag = list(range = TRUE))
}

\arguments{
  \item{x}{
A list with data tables as returned by \code{chop} or just one of such data tables (for \code{mf.flux} which will rarely be called directly by the user). Each table contains data for one chamber measurement. Usually there will be identifying columns, date, time and time sequence columns, and columns containing the concentrations of the greenhouse gases for each time step. At the minimum the columns specified in the argument \code{var.par} have to be there.
}
  \item{var.par}{
A named list specifying the variables and parameters that are used in the estimation process and variables that should be handed through the function so that they are easily available for further analysis. Some of the names are obligatory, others are optional. See details and examples.
}
  \item{subset}{
An optional vector specifying a subset of concentration measurements to be used in the estimation process.
}
  \item{asterisks}{
Logical. If TRUE p-values are given as asterisks.
}
  \item{loop}{
Logical or character. Set the mode how bulk estimations are run. \code{mf.flux} tries hard to find the linear part of the concentration development during chamber deployment. It does so by running many regressions. Therefore, when used in bulk mode with \code{fluxx} the estimation can be slow, especially when the data of many chamber deployments (contained in \code{x}) is computed. It may run faster when \code{loop} mode then. Per default (\code{"auto"}) the mode is changed automatically.
}
  \item{\dots}{
Further arguments passed to \code{mf.flux}.
}
  \item{time.unit}{
Single character giving the time unit. Will typically be seconds, thus default is \code{"S"}. Other options are \code{"M"} for minutes and \code{"H"} for hours. Be aware that fluxes will be incorrect when this setting does not reflect the time unit in your data.
}
  \item{all.through}{
Logical. When \code{TRUE}, all data columns in \code{x} other than the ones needed for flux calculation are also handed through the function so that they can be used in later steps of analysis easily. You may also specify one or a few columns that are handed through using \code{var.par}. The \code{all.through} setting overruns the \code{var.par} settings.
}
  \item{iv}{
Numeric. Sometimes there is no time information at all but the rows in \code{x} are just numbered consecutively. The correct temporal spread is calculated inside the function when the measurement interval is specified here. Defaults to 1 which expects times to be correctly given in \code{x}.
}
  \item{wndw}{
Numeric between 0 and 1. Relative width of a moving window in which the standard deviation of the concentrations is calculated to identify high frequency fluctuations. See details.
}
  \item{pdk}{
Numeric between 0 and 1. Minimum proportion of  data points to be kept. See details.
}
  \item{min.dp}{
Numeric. The minimum number of data points. Defaults to 20. If there are less rows the estimation is run anyway but a warning is issued and \code{min.dp} is automatically adjustet to n-1.
}
  \item{nrmse.lim}{
The maximum acceptable normalized root mean square error. Numeric value between 0 and 1. Defaults to 0.1. If the final best solution has a higher nrmse it is flagged accordingly.
}
  \item{r2.qual}{
Numeric between 0 and 1. Quality parameter for the model fit. The minimum acceptable r2 of the best fitted model. Defaults to 0.8. A measurement run is not defined invalid when the value is below the quality setting. Instead, a quality flag is reported.
}
  \item{range.lim}{
Numeric. The minimum range of the concentration measurements during one chamber placement. The acceptable range limit depends on the accuracy of the concentration measurements. When the range of the concentration measurements is smaller than the repeatability range of the measurement device one cannot tell real increase in concentration from random fluctuation. Therefore, if the range of the concentration measurements during one chamber placement is < \code{range.lim}, the range.flag is set to FALSE (0). See details.
}
  \item{out.unit}{
Output unit of the flux rate mass part. Character string. Defaults to "auto". In this case, the function tries to find a unit that ranges the output value between 0.01 and 10 which is fine for figures. Can be "ng", "mug", "mg", or "g". "mug" stands for "\eqn{\mu}g" because non-ascii characters are not allowed in functions. Beware of varying mass units in your output when running in auto mode.
}
  \item{elementar}{
When the fluxes are wanted as element values set \code{elementar = TRUE}. Defaults to FALSE.
}
  \item{hardflag}{
Named list that controls which of the quality flags are to be hard flagged (the value is changed according to the quality flag). Default is that the \code{range.lim} is hard flagged. So when the range of the actual concentration values for a chamber measurement is smaller than the set range limit a zero flux is given back. For changing a quality flag to a hard flag just provide the name of the quality parameter (without quotation marks) and set it to TRUE. Possible parameters that can be changed from/to hard/soft flag are \code{nrmse}, \code{range}, and \code{r2}. Further the number of measurements below ambient (\code{nomba}) can be hardflagged by setting \code{nomba = [0…nc]} with \code{nc} = number of concentration measurements for one chamber placement. Any flux measurement with a number of measurements below ambient >= \code{nomba} will be given back NA. 
}

}
\details{
The function is similar to \code{\link{flux}} but uses a different algorithm to identify the most linear part of the concentration development. First high frequency fluctations are omitted. Then all possible \code{pdk} * n : n consecutive concentration measurements are regressed against the corresponding times. The model with the highest r2 is chosen.

Probably the most important argument is \code{var.par}. It specifies the variables (by referring to the names of the data columns) from \code{x}  and  parameters (fixed values that are constant for all chamber placements) that are used for the flux estimations. The approach is similar to the one used in \code{\link{flux}} with the exception that there is no list item for the device quality flag. For obligatory \code{var.par} items see \code{\link{flux}} and examples.

In contrast to \code{\link{flux}} there is just one workhorse function doing the actual estimation (\code{mf.flux}) per data table. Especially when there are many data tables in \code{x} and/or many data points per data table it takes some time. Progress is shown in the console. Each dot represents one finalized data table.
}

\value{
\code{fluxx} returns a complex object of class \code{fluxxes} that is a 2 entry list. When the object is printed to the console only the second entry is displayed in a modified form that is meant to maximize information display with small footprint for easy inspection. A table is printed to the console with three columns per gas. The first contains the quality flags (e.g. 111.9). The order is: \code{r2.f}, \code{range.f},\code{nrmse.f}, \code{nomba.f}. The first three are considered more important, and if they are '1' everything is fine. For the one behind the dot its a bit different: It just gives the number of measurements below ambient.

The \code{\link[base]{data.frame}} with the estimated flux rates contains all data needed for further analysis. The columns represent the entries in \code{fluss} of the single chamber measurements (including quality flags, see below) plus naming information according to the settings in the \code{nmes} argument of \code{\link{chop}}. \code{\link{export}} provides a simple way to export the results.

The first entry is itself a list of lists and data tables. It is called \code{flux.res}. The only one first level entry in this list contains the information for one gas which is itself a list. In this list each first level entry contains the information for one chamber measurement. It is named according to the \code{nmes}-setting in \code{\link{chop}} and contains the elements \code{fluss} (which is itself a list with the elements given below), \code{mod}, \code{out} (a list with hand through data, list items according to columns in \code{x} that have been handed trough via \code{all.through} or \code{var.par}), and \code{inn} - a data.frame with the input data that were relevant for estimating the flux (the obligatory part of \code{var.par}).

The elements of \code{fluss}:

\item{ghg }{Character. The gas species for which the flux has been estimated.}
\item{flux }{Numeric. Calculated flux rate in mass unit per m2 and hour.}
\item{r2 }{r2 of the best fitted model that has been used for flux caclulation.}
\item{nrmse }{nrmse of the best fitted model that has been used for flux caclulation.}
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

\keyword{ univar }