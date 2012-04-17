\encoding{UTF-8}
\name{tt.nee}
\alias{tt.nee}
\alias{tt.flux}
\docType{data}
\title{
Medium frequency concentration data and fluxes from dynamic closed chamber measurements
}
\description{
The data comes from the Trebeltal / Northeastern Germany and has been recorded with flexible transparent and opaque dynamic closed chambers in 2011. 
}

\usage{
data(tt.nee)
data(tt.flux)
}

\format{
\code{tt.nee} is a \code{\link[base]{data.frame}} with 18 variables representing 14388 CO2 concentration measurements from 104 chamber placements.

\code{tt.flux} is a results table represnting fluxes estimated with \code{\link{fluxx}} from \code{tt.nee} with 28 columns and 104 rows (= number of chamber placements in \code{tt.nee}). Contains many variables from \code{tt.nee}. The additional ones start with item ???.
  
  \describe{
    \item{\code{date}}{Factor giving the date of field sampling, format is "\%Y-\%m-\%d".}
    \item{\code{time}}{Factor giving the time of measurement in the field, format is "\%H:\%M:\%S".}
    \item{\code{session}}{(Unique) Session number identifying one chamber placement. Integer.}
    \item{\code{record}}{Integer, running number of concentration measurement within one \code{session}. }  
    \item{\code{spot}}{Factor identifying the field measurement location.}
    \item{\code{PAR}}{Numeric. Photosynthetic photon flux density (PPFD).}
    \item{\code{t.cham}}{Numeric. Temperature logged inside chamber during concentration measurements}
    \item{\code{NEE}}{Numeric. CO2 concentration in chamber headspace.}
    \item{\code{t.air}}{Numeric. Air temperature outside chamber.}
    \item{\code{t.soil2}}{Numeric. Soil temperature at 2cm depth.}
    \item{\code{t.soil5}}{Numeric. Soil temperature at 5cm depth.}
    \item{\code{t.soil10}}{Numeric. Soil temperature at 10cm depth.}
    \item{\code{area}}{Numeric. Chamber area.}
    \item{\code{height}}{Numeric. Chamber height.}
    \item{\code{kind}}{Integer. Chamber kind. 1 = transparent chamber, 3 = transparent chamber with measurement before sun rise, 5 = opaque chamber}
    \item{\code{volume}}{Numeric. Chamber volume.}
    \item{\code{datetime}}{POSIXlt. Time stamp.}
    \item{\code{plot}}{Factor identifying the field plot (all TY1).}
    
    \item{\code{all}}{Factor. Combined unique identifier for chamber placement.}
    \item{\code{CO2.pv}}{Numeric. p.value of the fitted regression for the flux estimation.}
    \item{\code{CO2.r2.f}}{Logical numeric (0 | 1) giving the r2 flag. See \code{\link{fluxx}} for details.}
    \item{\code{CO2.range.f}}{Logical numeric (0 | 1) giving the range flag. See \code{\link{fluxx}} for details.}
    \item{\code{CO2.nrmse.f}}{Logical numeric (0 | 1) giving the nrmse flag. See \code{\link{fluxx}} for details.}
    \item{\code{CO2.ghg}}{Character. Greenhouse gas as submitted to \code{\link{fluxx}} via \code{var.par}.}
    \item{\code{CO2.unit}}{Character. Ouput unit of the flux as specified via \code{out.unit} in \code{\link{fluxx}}.}
    \item{\code{CO2.flux}}{Numeric. Flux}
    \item{\code{CO2.r2}}{Numeric. R2 of the model that has been used for flux estimation.}
    \item{\code{CO2.nrmse}}{Numeric. NRMSE of the model that has been used for flux estimation.}
    \item{\code{CO2.nomba.f}}{Numeric. Number of concentration measurements below ambient.}
    \item{\code{CO2.podpu}}{Numeric between 0 and 1. Propórtion of data points used.}
  }
}

\details{
\code{tt.nee} contains medium frequency (measured online) CO2 concentration data from 3 spots with \emph{Typha angustifolia} and includes data needed for modelling \code{\link{gpp}}/nne measured with transparent chamber and \code{\link{reco}} measured with opaque chamber. \code{tt.flux} contains fluxes estimated from \code{tt.nee} using \code{\link{fluxx}}. 
}

\source{
unpublished preliminary data
}

\examples{
## see examples at fluxx and gpp.
}

\keyword{datasets}
