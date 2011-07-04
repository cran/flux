\encoding{UTF-8}
\name{CO2.control}
\alias{CO2.control}
\title{
Simple function that controls how the CO2 concentration measurements influence the estimation of the flux rates for the other gases.
}

\description{
Just takes it arguments and hands them over to the argument \code{co2ntrol} in \code{\link{flux}}.
}

\usage{
CO2.control(CO2.columns = c("CO2ppm", "CO2Code"), trans.out = TRUE, range.lim = 30, min.allow = 3)
}

\arguments{
  \item{CO2.columns}{
A character vector of length 2 giving the names of the columns that contain the CO2 concentration measurements and the CO2 quality flag of the measurement device (see \code{\link{flux}} for details).
}
  \item{trans.out}{
Logical. Shall the CO2 outliers be transfered to \code{\link{flux.odae}} so that gas concentration measurements that caused CO2 outliers are not included in the optimation procedure anymore?
}
  \item{range.lim}{
Numeric value. The minimum ppm range for CO2 that is acceptable. For details see same argument at \code{\link{flux}}.
}
  \item{min.allow}{
Integer between 3 and the number of concentration measurements per chamber placement. For details see \code{min.allowed} at \code{\link{flux}}.
}
}

\details{
The idea behind \code{co2ntrol} in \code{\link{flux}} is that the CO2 concentration measurements might serve as a further check on the integrity of the chamber measuremnt in the field. For details see \code{\link{flux}}. With \code{trans.out = FALSE} it is possible to omit the transfer of CO2 outliers to the rest of the optimization procedure whilst maintaining the quality flag on leakage.
}

\value{
Returns a list with arguments that are used within \code{\link{flux}} to control the influence of CO2 concentration measurements on the flux estimation of other ghg. See arguments section for details.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\seealso{
\code{\link{flux}}, where you also find an example
}

\keyword{ univar }
