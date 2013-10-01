\encoding{UTF-8}
\name{flux}
\alias{flux}
\alias{flux.odae}
\alias{flux.conv}
\title{
Estimate gas flux rates using non-steady-state closed chamber data
}

\description{
\code{flux} is a convenience wrapper for the other two functions that should be suitable for most users. It can be used to estimate gas fluxes for all three commonly measured greenhouse gases (\eqn{CO_2}{CO2}, \eqn{CH_4}{CH4}, \eqn{N_2O}{N2O}) at once or separately.
}

\usage{
flux(x, var.par, co2ntrol = list(leak = TRUE, relay = FALSE), 
min.allowed = 3, max.nrmse = 0.1, nrmse.lim = 0.2, r2.qual = 0.8, 
range.lim = 30, out.unit = "auto", elementar = FALSE, 
hardflag = list(range = TRUE), asterisks = TRUE)

flux.odae(dat, var.par, min.allowed = 3, max.nrmse = 0.1, rl = NULL)

flux.conv(fl.dat, ghg = "CH4", r2.qual = 0.8, nrmse.lim = 0.2, 
out.unit = "auto", elementar = FALSE, hardflag = list(range = TRUE))
}

\arguments{
  \item{x}{
A list of data tables as returned by \code{chop}. Each table contains data for one chamber measurement. Required are at least the columns specified in the argument \code{var.par} (see also Example).
}
  \item{var.par}{
A named list specifying the variables and parameters that are used in the estimation process and variables that should be handed through the function so that they are easily available for further analysis. Some of the names are obligatory (\code{time}, \code{volume}, \code{area}, \code{t.air}, \code{p.air}, and two or more that specify the gas column and the gas quality column), others are optional. See details and examples.
}
  \item{ghg}{
Character string identifying the greenhouse gas for which concentration measurements are handled. Can be \code{"CH4"}, \code{"N2O"} or \code{"CO2"}. Defaults to \code{"CH4"}. This argument determines the molar weight that is used in the flux calculation with \code{\link{gflux}}, the input unit (ppb or ppm), and the ambient concentration of the gas that is added to the diagnostic plots plotted with \code{\link{plot.fluss}} and used for the determination of the number of measurements below ambient (\code{nomba}).
}
  \item{co2ntrol}{
Options for estimating fluxes with \eqn{CO_2}{CO2} control. In this case outliers and the slope of the \eqn{CO_2}{CO2} concentration measurements influence the estimated flux rate of the other greenhouse gases (\eqn{N_2O}{N2O} and \eqn{CH_4}{CH4}; see Details). By default, \eqn{CO_2}{CO2} correction is switched off.
}
  \item{min.allowed}{
Integer giving the minimum number of concentration measurements allowed during the estimation of one single flux. Can be any number between 3 and the number of concentration measurements during one chamber placement.
}
  \item{max.nrmse}{
Numeric giving the maximum acceptable normalized root mean square error for configurations with higher numbers of concentration measurements than specified in \code{min.allowed}. Numeric value between 0 and 1. Defaults to 0.1. Above that value lower numbers of concentration measurements down to \code{min.allowed} are considered. See details. In \code{flux} \code{max.nrmse} can also be given as a named list with three items giving the maximum acceptable nrmse per gas. See examples.
}
  \item{nrmse.lim}{
Numeric between 0 and 1 (defaults to 0.2) giving the main quality parameter for the model fit, the maximum acceptable normalized root mean square error. If the best fit for one chamber measurement exceeds this value, the function reports \code{FALSE} in the \code{nrmse.f} quality flag. See details and value. In \code{flux} \code{nrmse.lim} can also be given as a named list with three items giving the \code{nrmse.lim} per gas. See examples.
}
  \item{range.lim}{
The minimum detectable range of the concentration measurements during one chamber placement. Has to be either a single numerical value, a numeric vector with the same length as \code{x} giving different range limits for each chamber placement (for instance obtained by \code{\link{flux.calib}}) or a character string naming the column in \code{x} that contains range limit data. If this column is named "rl" (the default when the range limits are attached to the data by \code{\link{flux.calib}})  \code{flux} automatically detects it. Note, that setting \code{range.lim != NULL} overrides the auto detected range limits.

The acceptable range limit depends on the accuracy of the concentration measurements. When the range of the concentration measurements is smaller than the repeatability range of the measurement device (e.g., a gas chromatograph) one cannot tell real increase in concentration from random fluctuation. Therefore, if the range of the concentration measurements during one chamber placement is < \code{range.lim}, the range.flag is set to FALSE (0). See details. In \code{flux} \code{range.lim} can also be given as a named list with three items giving the maximum acceptable range limits per gas. See examples.
}
  \item{r2.qual}{
  Numeric giving the limit of minimum acceptable r2 as an alternative quality parameter describing the model fit. Can be between 0 and 1 (0.8 by default). If a model r2 is below the setting the r2 quality flag is reported FALSE (0). In \code{flux} \code{r2.qual} can also be given as a named list with three items giving the acceptable \code{r2.qual}ities per gas. See examples.
}
  \item{out.unit}{
Character string determining the output unit of the flux rate mass part. Character string. The default "auto" tries to find a unit that ranges the output value between 0.01 and 10. Possible output units are "ng", "mug", "mg", or "g". "mug" stands for "\eqn{\mu}g" because non-ascii characters are not allowed in functions.
}
  \item{elementar}{
When the fluxes are wanted as element values set \code{elementar = TRUE}. Defaults to FALSE.
}  	
  \item{hardflag}{
Named list that controls which of the quality flags are to be hard flagged (the value is changed according to the quality flag). \code{range.lim} is hard flagged by default. So when the range of the actual concentration values for a chamber measurement is smaller than the set range limit a zero flux is returned. When a flux estimation does not meet the quality requirements of any other hard-flagged quality parameter \code{flux} returns NA. For changing a quality parameter to a hard flag just provide its name (without quotation marks) and set it to TRUE. Possible parameters to chose as hard or soft flag are \code{nrmse}, \code{range}, and \code{r2}. Further the number of measurements below ambient (\code{nomba}) can be hardflagged by setting \code{nomba = [0…nc]} with \code{nc} = number of concentration measurements during one chamber placement. See examples.
}
  \item{asterisks}{
Logical. If \code{TRUE}, p-values are given as asterisks.
}
  \item{dat}{
One data table for one chamber placement. See \code{x} and \code{var.par} for details.
}
  \item{rl}{
Specifies \code{range.lim} in the low level function \code{flux.odae}. As with \code{range.lim} several options are allowed. Defaults to NULL. In this case, the function looks for a column \code{rl} in \code{dat}. If it can't find a column \code{rl}, the value is set to 0 and a warning is returned; if it does exist \code{dat$rl} is always taken. If \code{rl} is a character string the function looks for a column of that name in \code{dat}, if \code{rl} is a numeric value, this value is taken as the range limit. See examples.
}
  \item{fl.dat}{
An object with the same structure as returned by \code{\link{flux.odae}}. See details and value.
}
}

\details{
Typically it will be most convenient to use \code{flux} on objects returned by \code{\link{chop}} (i.e. on lists of data tables that contain all necessary data per chamber measurement including supporting information). \code{flux} simply wraps \code{flux.odae} and \code{flux.conv} applied on lists of chamber measurement data tables into one function. Thus, the data of a one day field campaign or a year of chamber measurements can easily be handled by simply running two functions (\code{chop} and \code{flux}) consecutively to estimate ghg fluxes for all three common ghg gases. See example.

Probably the most important argument is \code{var.par}. It specifies the variables (by referring to the names of the data columns) from \code{x}  and  parameters (fixed values that are constant for all chamber placements) that are used for the flux estimations. For simple handling it is expecting a named list. 

For \code{flux} the obligatory list items are: One item that refers to a column in \code{x} containing ghg concentrations (see next paragraph for details); \code{time} – chamber closing time in minutes; \code{volume} – chamber volume during placement (in cbm); \code{area} – chamber area (in sqm); \code{t.air} – air temperature inside chamber during concentration measurements (in °C); \code{p.air} – air pressure during concentration measurements (in Pa). 

The list items that are used to specify the ghg for which flux estimation is carried out have to be specified by using the named list items \code{CH4} – \eqn{CH_4}{CH4} concentrations; \code{CH4.gcq} – \eqn{CH_4}{CH4} GC quality flag; \code{CO2} – \eqn{CO_2}{CO2} concentrations; \code{CO2.gcq} – \eqn{CO_2}{CO2} GC quality flag; \code{N2O} – \eqn{N_2O}{N2O} concentrations; \code{N2O.gcq} – \eqn{N_2O}{N2O} GC quality flag. Fluxes are estimated for all ghg for which concentration data are given. Thus, at least one ghg should be specified. GC quality flags are optional. If you don't provide a reference to a column in \code{x} the function assumes that all GC measurements were OK.

All these list items can either be given as a variable (name of a column in \code{x}) or as a fixed parameter (a numeric value). This makes no sense for the \code{ghg}s and \code{time}, but in many cases chamber \code{volume} and \code{area} will be constant across measurements. Another likely candidate for a fixed parameter is \code{p.air} because air pressure is often not logged during chamber measurements. All additional list items should be of type ‘variable‘ and refer to further columns in \code{x} if you want those data handed through the function and be part of the result tables (for having all data in one place for further analyses). You are free to choose appropriate names. Fixed parameters will not be relayed.

If the flux estimation is carried out in two steps it will typically be carried out on a list structure as returned by \code{\link{chop}}. Therefore, it is used within a \code{\link{lapply}} call. For details see examples. However, the functions \code{flux.odae} and \code{flux.conv} are  designed to be carried out on single data tables (\code{\link[base]{data.frame}}) per chamber measurement. 

First \code{flux.odae} is run. It simply tries to find the best model fit for the series of concentration measurements that are given in \code{dat}. This \code{\link[base]{data.frame}} has to consist of five columns that give (in that order): gas concentration, closing time of the chamber in minutes, gas concentration quality flag, chamber volume, temperature within the chamber headspace during measurements (may change during chamber placement). See example data. 

At the moment the optimization bases on linear regression. All possible models with n (= total number of concentration measurements per chamber placement) to \code{min.allowed} number of concentration measurements are fitted and the best fit is evaluated in a stepwise procedure. The normalized root mean square error is used as the quality criterion for the outlier detection and elimination procedure. All model fits with a nrmse <= \code{max.nrmse} are extracted and ranked according to the number of concentration measurements (decreasing) and to the nrmse (increasing). The first ranked model is stored along with the original data table and some other information. Therefore a model with e.g. a nrmse of 0.081 constructed from 5 concentration measurements wins against a model with a nrmse of 0.07 with only 4 concentration measurements. This reflects the idea that models with nrmse <= \code{max.nrmse} already represent a sufficient fit and do not have "outliers" that must be eliminated. 

In case no model has a nrmse <= \code{max.nrmse}, the models are simply ranked according to their nrmse and the model with the lowest nrmse wins and is stored. In that way outliers are detected and exluded. \code{flux.odae} returns a complex object that contains most of the necessary information for the \code{flux.conv} function and also carries information that is later needed for the plot functions (\code{plot.flux} and \code{plot.fluss}). 

The flux calculation is then carried out with the function \code{flux.conv}. It takes the object returned by \code{flux.odae} and additional information (chamber area, gas species, several quality settings and in- as well as output units) and calculates the flux rates. Further several quality checks (r2 check, range check, nrmse check, nomba check; for details see Value) are carried out and quality flags are reported along with the fluxes in the output. It is best when all quality flags are returned \code{TRUE}. Depending on the application quality requirements might vary. Therefore, per default the function reports soft quality flags (despite for range). However, this can be changed via \code{hardflag}.

The idea behind \code{co2ntrol} in \code{\link{flux}} is that the \eqn{CO_2}{CO2} concentration measurements might serve as a further check on the integrity of the chamber measuremnt in the field. When \code{co2ntrol} is set, the function first carries out an outlier procedure on the \eqn{CO_2}{CO2} concentration data (the respective columns have to be in \code{x} of course). Further, the slope of the \eqn{CO_2}{CO2} concentration change over time is checked. When it is negative, chamber leakage is assumed and a respective quality flag (\code{leak.flag}) is reported FALSE. The \code{leak.flag} cannot be hard flagged.
}

\value{
\code{flux} returns a complex object of class \code{fluxes} that is a 3 entry list. When the object is printed to the console only the second entry is displayed in a modified form that is meant to maximize information display with small footprint for easy inspection. A table is printed to the console with three columns per gas. The first contains the quality flags (e.g. "111.02"). The order is: \code{nrmse.f}, \code{r2.f}, \code{range.f}, \code{nomba.f}, \code{leak.f}. The first three are considered more important, and if they are '1' everything is fine. The first flag behind the full stop just gives the number of measurements below ambient, while the second is '2' when \code{co2ntrol} was switched off, '0' when leaking occurred, and '1' when no leaking occurred.

The \code{\link[base]{data.frame}} with the estimated flux rates contains all data needed for further analysis. The columns represent the entries in \code{fluss} of the single chamber measurements (including quality flags, see below) plus naming information according to the settings in the \code{nmes} argument of \code{\link{chop}}. \code{\link{export}} provides a simple way to export the results. The first entry is itself a list of lists and data tables. It is called \code{flux.res} and is comprised of objects that are returned by \code{flux.conv} per ghg. Each first level entry in these lists contains the information for one chamber measurement. It is named according to the \code{nmes}-setting in \code{\link{chop}} and contains the elements \code{fluss} (which is itself a list with the elements given below), \code{fl.dat} (equals the object returned by \code{flux.odae}; see below), and \code{unit} which provides information on the output mass unit of the flux rate that is handed over to the function \code{\link{plot.fluss}} and to the table output.

The elements of \code{fluss}:

\item{ghg }{Character. The gas species for which the flux has been estimated.}
\item{flux }{Numeric. Calculated flux rate in mass unit per m2 and hour.}
\item{r2 }{The \eqn{R^2}{R2} of the best fitted model that has been used for flux caclulation.}
\item{nrmse }{The NRMSE of the best fitted model that has been used for flux caclulation.}
\item{r2.f }{Logical. \eqn{R^2}{R2} quality flag telling whether the \eqn{R^2}{R2} quality setting given in \code{r2.qual} is fulfilled.}
\item{range.f }{Logical. Range quality flag telling whether the range of the concentration measurements exceeded the quality range of the measurement device that has been specified in \code{range.lim}.}
\item{nrmse.f }{Logical. NRMSE quality flag telling whether the NRMSE quality setting given in \code{nrmse.lim} is fulfilled (i.e. if the NRMSE of the best model <= \code{nrmse.lim}).}
\item{nomba.f }{Integer. Reports the \strong{n}umber \strong{o}f \strong{m}easurements \strong{b}elow \strong{a}mbient. The ambient concentrations are set to be 392.6 ppm (\eqn{CO_2}{CO2}), 1874 ppb (\eqn{CH_4}{CH4}),  and 324 ppb (\eqn{N_2O}{N2O}) (taken from Mace Head Ireland (\eqn{N_2O}{N2O}, \eqn{CH_4}{CH4}) and global average (\eqn{CO_2}{CO2}) obtained from http://cdiac.ornl.gov/pns/current_ghg.html as of August 16th, 2013).}
\item{leak.f }{Logical. When \code{co2ntrol} was applied with \code{leak = TRUE}, possible chamber leakage as represented by decreasing \eqn{CO_2}{CO2} concentrations over time is shown by a \code{FALSE (0)}.}

The elements of \code{fl.dat} that is also the object returned by \code{flux.odae} are:

\item{lm4flux }{Complex object. The best fitting model as reported by \code{\link{lm}}. It builds the basis for the calculation of the flux rate via \code{flux.conv}.}
\item{row.select }{Integer vector giving the indices of the rows of the data table that have been used to construct the best fitting model. This information is later used in the plotting functions \code{\link{plot.flux}} and \code{\link{plot.fluss}}. }
\item{orig.dat }{\code{\link[base]{data.frame}} with the original data provided according to arguments \code{x} and \code{columns}. }
\item{out.dat }{Data to be handed through. Per default \code{area} and \code{volume} of the chamber are relayed but these values are not part of the table output whereas all additionally relayed data are part of the table output. }

}

\references{
Nakano T (2004) A comparison of regression methods for estimating 
soil-atmosphere diffusion gas fluxes by a closed-chamber technique. 
Soil Biology and Biochemistry 36: 107-113.

Forbrich I, Kutzbach L, Hormann A, Wilmking M (2010) A comparison of 
linear and exponential regression for estimating diffusive CH4 fluxes by 
closed-chambers in peatlands. Soil Biology and Biochemistry 42: 507-515.
}

\author{ Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>, Franziska Koebsch <franziska.koebsch@uni-rostock.de>, Ulrike Hagemann <ulrike.hagemann@zalf.de>, Anke Günther <anke.guenther@uni-rostock.de>
}

\seealso{
\code{\link{chop}}, \code{\link{flux.calib}}, \code{\link{gflux}}, \code{\link{plot.fluss}}
}

\examples{
## load example data
data(tt.pre)

## extract field concentration measurements
gcd <- tt.pre[tt.pre$sampletype_a=="P",]

## partition the data into data tables per chamber measurement
# then do the partitioning
gcd.parts <- chop(gcd, factors = c("date", "spot", "veg"), 
nmes = c("date", "veg", "spot"))

## calculate flux rates for methane
# first define a global CH4 range limit
CH4.lim <- 30
# do the flux rate estimation (it will often be best to define
# var.par separately, note that p.air is given as a parameter)
vp.CH4 <- list(CH4 = "CH4ppb", time = "time_min", CH4.gcq = "CH4Code", 
volume = "cham_vol", t.air = "temp_dC", area = "cham_area", p.air = 101325)
flux.CH4 <- flux(gcd.parts, var.par = vp.CH4)
# look at the results table
flux.CH4

# extracting range limits from the calibration gas measurements
# and attaching them to gcd.parts. first get the calibration gas
# measurements from tt.pre (changing the date because it is in 
# a strange format and has to be the same as the dates in gcd.parts)
cgm <- tt.pre[tt.pre$sampletype_a=="E",c("date_gc", "CH4ppb", "CH4Code", 
"CO2ppm", "CO2Code", "N2Oppb", "N2OCode")]
names(cgm)[1] <- "date"
cgm$date <- "2011-03-16"
# now we can do the flux.calib
gcd.parts.cal <- flux.calib(gcd.parts, columns = c("date", "CH4ppb"), 
calib = cgm, format="\%Y-\%m-\%d", window=48, buffer=1100, attach=TRUE)
# do the flux rate estimation (we use the same var.par as before)
flux.CH4 <- flux(gcd.parts.cal, var.par=vp.CH4, co2ntrol = NULL, 
range.lim=NULL)
# look at the results table
flux.CH4
# export the results to the working directory
wd <- getwd()
export(flux.CH4, file=paste(wd, "/flux.CH4.txt", sep=""))

## plot the concentration-change-with-time-plots as kind of diagnostic
plot(flux.CH4, dims = c(3,6))

## do the flux rate estimation whilst using CO2 concentrations to
## control for possible chamber leakage
flux.CH4.b <- flux(gcd.parts, var.par=vp.CH4)
# look at the results table
flux.CH4.b
# plot the concentration-change-with-time-plots as kind of diagnostic
plot(flux.CH4.b, dims = c(3,6))

## do the flux rate estimation whilst using CO2 concentrations to
## control for outliers and possible chamber leakage
flux.CH4.c <- flux(gcd.parts, var.par=vp.CH4, co2ntrol = list(leak = TRUE, 
relay = FALSE))
# look at the results table
flux.CH4.c
# plot the concentration-change-with-time-plots as kind of diagnostic
plot(flux.CH4.c, dims = c(3,6))

}

\keyword{ univar }
