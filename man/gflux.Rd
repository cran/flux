\encoding{UTF-8}
\name{gflux}
\alias{gflux}
\title{
Calculate gas flux rate from two concentrations
}
\description{
Calculate gas flux rate from two concentrations using the ideal gas law to obtain a mass flow from an area per time. Therefore, besides the two concentrations ct and c0 the temperature within and the volume of the closed chamber are needed. For areal reference the area from which the gases are emitted has to be given. Without any further unit transformation the input unit directly gives the output unit: When concentration is coming in ppm the calculated flux rate is in \eqn{\mu}g/m2*h, when concentration is in ppb the flux rate will be in ng/m2*h
}
\usage{
gflux(ct, c0 = NULL, T, V, A = 1, M = 44, t = 1/60, p = 101325)
}

\arguments{
  \item{ct}{
Concentration of the gas at time t. When the function is used internally the concentration is automatically derived from a regression model. May also be the concentration change in time dc/dt. In this case \code{c0} must not be specified.
}
  \item{c0}{
Concentration of the gas at time 0. When the function is used internally the concentration is automatically derived from a regression model.
}
  \item{T}{
Temperature within the chamber during the measurement in \eqn{^{\circ}}C (it is converted automatically to Kelvin). When it is changing during the measurement typically the average temperature is used. However, if there is too much change in temperature during the chamber closing time (more than 5 \eqn{^{\circ}}K) the ideal gas law might not longer be appropriate.
}
  \item{V}{
Chamber headspace volume. Because concentrations are typically small volume matters and should therefore be determined as exactly as possible.
}
  \item{A}{
Area covered by the measurement chamber. Defaults to 1. For dimensionless sampling just leave the default. As with the volume it matters a lot for the end value, therefore it should be determined as exactly as possible.
}
  \item{M}{
Molar weight of the gas for wich concentration data is given. Defaults to 44 because two of the three most commonly considered greenhouse gases share this molar weight (\eqn{N_2O}{N2O} and \eqn{CO_2}{CO2}). When calculating the flux rate from methane concentrations change accordingly (16 g/mol).
}
  \item{t}{
Chamber closing time or more exactly the time span between the measurement of c0 and ct. When derived from a regression model this might not be the whole chamber closing time.
}
  \item{p}{
The air pressure at earth surface during measurements. Default is given by the standard value at sea level of 101325 Pa. Should be OK for most lowland measurements. However, if the measurements took place on higher altitudes, it might be reasonable to adapt the air pressure value to local conditions.
}
}

\details{
Typically the function will not be called separately. However, for checking single chamber measurements or for testing purposes it might be useful to have this as a separate function.

The flux rate is calculated using

flux.rate = ((ct-c0) * V * M * p) / (t * R * (T + 273.15) * A)

The gas constant R is used with its standard value 8.314 Pa/K*mol.
}

\value{
Returns one numeric value that represents the flux rate in mass unit (depending on input concentration) per \eqn{m^2}{sqm} and hour.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>, 
Franziska Koebsch <franziska.koebsch@uni-rostock.de>
}

\seealso{
\code{\link{flux}}, \code{\link{flux.conv}}
}

\examples{
## see examples for function flux
}

\keyword{ univar }
