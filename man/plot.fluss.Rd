\encoding{UTF-8}
\name{plot.fluss}
\alias{plot.fluss}
\alias{plot.flux}
\alias{plot.fluxes}
\alias{plot.fluxx}
\alias{plot.fluxxes}
\title{
Functions for bulk plotting of concentration change with time as kind of 
diagnostic plots for flux rate calculations
}
\description{
Bulk plotting of concentration change with time adding color and symboling  
for acting as a diagnostic plot for the flux rate estimation functions 
(\code{\link{flux}}, \code{\link{flux.odae}}, \code{\link{flux.conv}}) in 
this package.
}

\usage{
\method{plot}{fluss}(x, subs, dims, folder = getwd(), xlims = NULL, ...)

\method{plot}{flux}(x, zero.line, note = "", margin = 0.2, xlims = NULL, ...)

\method{plot}{fluxes}(x, dims, ghg = "all", subs = NULL, folder = getwd(), 
xlims = NULL, ask = TRUE, ...)

\method{plot}{fluxx}(x, ...)

\method{plot}{fluxxes}(x, dims, subs = NULL, folder = getwd(), ask = TRUE, ...)
}

\arguments{
  \item{x}{
Object of class \code{fluss} that is returned by \code{\link{flux}} or object of class \code{flux} that is returned by \code{\link{flux.conv}}. In case of the latter, the function has to be applied to a list of flux extimation results via \code{\link[base]{lapply}} or the like.
}
  \item{subs}{
Single character value or character value specifying the factors that shall be used for subsetting the plots into plates (that are stored as pdf files to a folder specified in \code{folder}). Must be names of columns of the original data that have been used to partition the data into chunks and that are part of the naming of the data chunks (see \code{\link{flux}} for details.). When there are only few chamber measurements that shall be plotted to the screen set \code{subs} to NULL. This is default behaviour for \code{plot.fluxes}
}
  \item{dims}{
Integer vector with two elements that specify the mfrow setting (see \code{\link[graphics]{par}} for details) during the plotting of the single plates into pdf files. For all single plots to fit on the plate the product of the two entries has to be equal or higher the number of chamber measurements that are in the data in each partition according to \code{subs}.
}
  \item{folder}{
Character string giving the path to the folder were the files have to be stored. The names of the pdf files are generated automatically.
}
  \item{xlims}{
Two entry numeric vector specifying the x-axes limits for all plots. Defaults to NULL in which case it is derived from the data itself. The y-axes limits are always set according to the range of the concentration data \eqn{\pm} 20% (this can be changed via \code{margin} if necessary). 
}
  \item{\dots}{
further arguments passed through to \code{plot.flux} (see below) or to \code{\link[graphics]{plot.default}}.
}
  \item{zero.line}{
The y-axes position of a horizontal line that reflects the ambient concentration of the plotted gas species. When using \code{plot.fluss} this is determined automatically from \code{x}.
}
  \item{note}{
A note that shall appear in the plots. Typically not a fixed value but a value that changes from plot to plot. See example.
}
  \item{margin}{
Numeric between 0 and 1. Specifies the empty space within the diagnostic plots on the y-axis. 
}
  \item{ghg}{
Character value or an up to three entry vector specifying which ghg should be plotted. Note that only ghg fluxes that were estimated can be plotted.
}
  \item{ask}{
Logical; if TRUE, the user is asked before starting to plot the concentration data for the next ghg, see \code{\link[graphics]{par}}(ask=.) and examples.
}
}

\details{
Typically \code{plot.fluss} will be used. However, for lower level plotting the function \code{plot.flux} that also does the plotting within \code{plot.fluss} is provided as a separate function.
}

\value{
The function is invoked for its side effects and does not return anything.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\seealso{
\code{\link{chop}}, \code{\link{flux}}, \code{\link{flux.odae}}, \code{\link{flux.conv}}
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
# first define CH4 range limit (alternatively use flux.calib)
CH4.lim <- 30
# do the flux rate estimation
vp.CH4 <- list(CH4 = "CH4ppb", time = "time_min", CH4.gcq = "CH4Code", 
volume = "cham_vol", t.air = "temp_dC", area = "cham_area", p.air = 101325)
flux.CH4 <- flux(gcd.parts, var.par = vp.CH4, co2ntrol = NULL, 
range.lim=CH4.lim)

## look at the results table
flux.CH4

## plot the concentration-change-with-time-plots as kind of diagnostic
plot(flux.CH4, dims = c(3,6))
}

\keyword{ hplot }