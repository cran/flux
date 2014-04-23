\encoding{UTF-8}
\name{flux-package}
\alias{flux-package}
\docType{package}
\title{
Flux rate estimation with dynamic closed chamber data
}
\description{
Several functions for the estimation of greenhouse gas (GHG) flux rates using closed chamber concentration measurements. The package follows a modular concept: Fluxes can be calculated in just two simple steps or in several steps if more control is wanted. Functions for further analyses (GPP and Reco model fitting and prediction for budgets including error terms) are also available.
}
\details{
\tabular{ll}{
Package: \tab flux\cr
Type: \tab Package\cr
Version: \tab 0.3-0\cr
Date: \tab 2014-04-23\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Obtain flux rates from many chamber measurements within minutes. After 
preparing the read in data (Field or device measured concentration data on 
the three most prominent greenhouse gases) with \code{\link{chop}} just 
run \code{\link{flux}} or \code{\link{fluxx}} (for medium frequency data)
on the result returned by \code{\link{chop}} and get flux rates in an 
easy to interpret table including quality flags. 
Plot diagnostic plots as pdf per factor level to a folder or simply to 
the screen. Use \code{\link{gpp}} to model GPP and \code{\link{reco}} to 
model ecosystem respiration or use \code{\link{gpp.bulk}} to bulk model
GPP and \code{\link{reco.bulk}} to bulk model \eqn{R_{eco}}{R_eco} and 
use the resulting objects with \code{\link{budget.gpp}} and 
\code{\link{budget.reco}}, respectively, to predict fluxes using 
continuously logged data. Use \code{\link{budget.ie}} to estimate the 
uncertainty associated with the interpolation between models. 
Several helper functions for ghg analysis are also provided.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>, Franziska Koebsch <franziska.koebsch@uni-rostock.de>, Ulrike Hagemann <ulrike.hagemann@zalf.de>, Anke Günther <anke.guenther@uni-rostock.de>

Maintainer: Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\references{
Nakano T (2004) A comparison of regression methods for estimating soil-
atmosphere diffusion gas fluxes by a closed-chamber technique. Soil Biology 
and Biochemistry 36: 107-113.

Forbrich I, Kutzbach L, Hormann A, Wilmking M (2010) A comparison of linear 
and exponential regression for estimating diffusive CH4 fluxes by closed-
chambers in peatlands. Soil Biology and Biochemistry 42: 507-515.

Beetz S, Liebersbach H, Glatzel S, Jurasinski G, Buczko U, Hoper H (2013) Effects of land use intensity on the full greenhouse gas balance in an Atlantic peat bog. Biogeosciences 10:1067-1082.

Koebsch F, Glatzel S, Jurasinski G (2013) Vegetation controls methane emissions in a coastal brackish fen. Wetlands Ecology and Management 21:323–337.
}

\keyword{ package }

\seealso{
\code{\link[HMR]{HMR}} for a different approach to flux rate estimation from chamber data.
}
