\encoding{UTF-8}
\name{modjust}
\alias{modjust}
\title{
Adjust Reco models
}
\description{
The function allows to adjust fitted \eqn{R_{eco}}{R_eco} models by eliminating the maximum \eqn{R_{eco}}{R_eco} flux as long as the p.value of the linear model of the residuals regressed against original fluxes is above a given threshold. In addition models with parameters that went astray may be skipped. The default is that \eqn{R_{eco}}{R_eco} models with t1 > 20 are omitted.
}

\usage{
modjust(models, alpha = 0.1, minimum = 0.8, prmtrs = list(t1 = 20), ...)
}

\arguments{
  \item{models}{
Object of class "\code{breco}".
}
  \item{alpha}{
Alpha level against which the p.value of the linear model of the residuals against original fluxes shall be tested.
}
  \item{minimum}{
The minimum proportion of data points that should be kept. The optimisation runs in a \code{while} loop until the p.value is below \code{alpha}. It may happen - especially when the number of data points was already low from beginning - that many data points are skipped before a solution is reached. This is prevented by this argument which acts as a brute force to the process and stops it. 
}
  \item{prmtrs}{
List object that allows to skip models according to thresholds set for coefficients of the fitted regression models. The list has to be set up according to the actual method used in \code{\link{reco}} an the names refer to the names of the corresponding coefficients. The default is that \eqn{R_{eco}}{R_eco} models with t1 > 20 are omitted.
}
  \item{\dots}{
Arguments passed through to \code{\link{reco}} which is used to fit the models again based on the adjusted data.
}
}

\details{
When fitting \eqn{R_{eco}}{R_eco} models based on one or few measurement campaigns in the field it may happen that outliers in the extremes of the temperature gradient have a very high influence on the fit. Although the model could be fit in the first place this often leads to unrealistic predicted fluxes. The adjustment via \code{modjust} leads to better overall performance and reliability of the bulk modelling.
}

\value{
Returns a "\code{breco}" object with the possibly adjusted models. All returned models gain a list entry within the \code{mod} object (see \code{\link{reco}} and \code{\link{reco.bulk}}) named \code{n.out.adj} giving the number of omitted data points. Fall.back models (see \code{\link{reco.bulk}}) in \code{models} are left untouched.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de},

based on ideas by Sascha Beetz, \email{sascha.beetz@uni-rostock.de}
}

\seealso{
\code{\link{reco}}, \code{\link{reco.bulk}}
}

\examples{
## See axamples at reco.bulk
}

\keyword{ manip }
