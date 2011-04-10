\encoding{UTF-8}
\name{flux-package}
\alias{flux-package}
\docType{package}
\title{
Flux rate estimation with dynamic closed chamber data
}
\description{
Several functions for the estimation of greenhouse gas flux rates using dynamic closed chamber concentration measurements. The package follows a modular concept: Fluxes can be calculated in just two simple steps or in several steps if more control is wanted.
}
\details{
\tabular{ll}{
Package: \tab flux\cr
Type: \tab Package\cr
Version: \tab 0.1-2\cr
Date: \tab 2011-04-10\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Obtain flux rates from many chamber measurements within minutes. After 
preparing the read in data (Field or device measured concentration data on 
the three most prominent greenhouse gases) with \code{\link{conz.prep}} just 
run \code{\link{flux}} on the result returned by \code{\link{conz.prep}} and 
get flux rates in an easy to interpret table including quality flags. Plot 
diagnostic plots as pdf per factor level to a folder or simply to the 
screen.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>, Franziska Koebsch <franziska.koebsch@uni-rostock.de>

Maintainer: Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\references{
Nakano T (2004) A comparison of regression methods for estimating soil-
atmosphere diffusion gas fluxes by a closed-chamber technique. Soil Biology 
and Biochemistry 36: 107-113.

Forbrich I, Kutzbach L, Hormann A, Wilmking M (2010) A comparison of linear 
and exponential regression for estimating diffusive CH4 fluxes by closed-
chambers in peatlands. Soil Biology and Biochemistry 42: 507-515.
}

\keyword{ package }

\seealso{
\code{\link[HMR]{HMR}} for a different approach to flux rate estimation from chamber data (aditionally it is meant for static chamber data).
}
