\encoding{UTF-8}
\name{budget.ie}
\alias{budget.ie}
\title{
Estimate interpolation errors for GPP and Reco budgets
}
\description{
Use several MC methods to estimate the uncertainty associated with the interpolation of fluxes between models when preparing time series data for budgeting \eqn{R_{eco}}{R_eco}, GPP (and finally NEE) with \code{\link{budget.reco}} and \code{\link{budget.gpp}}. 
}

\usage{
budget.ie(bdgt, method = "leave out", lo = 2, it = 100)
}

\arguments{
  \item{bdgt}{
An object deriving from running \code{\link{budget.reco}} or \code{\link{budget.gpp}}. The budget should cover a reasonable time span, e.g. a year or even better, a calender year and must have been derived with \code{return.models} = \code{TRUE} and with a \code{set.back} that only defines start and end dates. See example at \code{\link{budget.reco}}.
}
  \item{method}{
Specify how interpolation error should be estimated. Available methods include \code{"leave out"}, \code{"bootstrap"}, \code{"sorted bootstrap"}, \code{"constrained bootstrap"}, \code{"jackknife"}, \code{"jack-validate"}. True bootstrap is only effective when \code{sort.x} = \code{FALSE}. See \code{\link{auc.mc}} for details.
}
  \item{lo}{
When estimating interpolation error with \code{"leave out"} or \code{"jack-validate"}, how many data points should be left out randomly? Defaults to 2. See \code{method} and details.
}
  \item{it}{
How many iterations should be run? Defaults to 100. Not effective for methods \code{"jackknife"}, and \code{"jack-validate"}.
}
}

\details{
ATTENTION: It takes a while. How long one budget run takes depends on the length of the \code{bdgt} but typically takes about 5 seconds. So if you run with defaults (\code{it} = 100) it may take some minutes. Progress is shown in the console with numbers representing the runs separated by colons.

The approaches are quite similar to the ones in \code{\link{auc.mc}}. However, the function randomly samples from a list of models and then runs the complete budgeting via \code{\link{budget.reco}} or \code{\link{budget.gpp}}. This is done either \code{it} times or as often as needed to get all combinations that are possible (for methods \code{"jackknife"}, and \code{"jack-validate"}).
}

\value{
A vector of budgets. 
}

\references{
Beetz S, Liebersbach H, Glatzel S, Jurasinski G, Buczko U, Hoper H (2013) Effects of land use intensity on the full greenhouse gas balance in an Atlantic peat bog. Biogeosciences 10:1067-1082
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de},
}

\seealso{
\code{\link{reco.bulk}}, \code{\link{gpp.bulk}}, \code{\link{modjust}}, \code{\link{budget.reco}}, \code{\link{budget.gpp}}
}

\examples{
## See examples at budget.reco

}

\keyword{ univar }
