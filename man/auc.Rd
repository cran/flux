\encoding{UTF-8}
\name{auc}
\alias{auc}
\alias{auc.mc}
\title{
Calculate the area under a line(curve).
}
\description{
Calculates the \strong{a}rea \strong{u}nder a \strong{c}urve (integral) following the trapezoid rule. With \code{auc.mc} several Monte Carlo methods can be applied to obtain error terms for estimating the interpolation error for the integration.
}

\usage{
auc(x, y, thresh = NULL, dens = 100, sort.x = TRUE)

auc.mc(x, y, method = "leave out", lo = 2, it = 100, ...)
}

\arguments{
  \item{x }{ Numerical vector giving the x cordinates of the points of the line (curve). }
  \item{y }{ Numerical vector giving the y cordinates of the points of the line (curve). One can calculate the integral of a fitted line through giving a vector to \code{x} that spans \code{xlim} with small intervals and predicting the y coordinates with \code{\link[stats]{predict}} and that \code{x}-vector as \code{newdata}. See example. }
  \item{thresh }{ Threshold below which area is not calculated. When below threshold data represent proper data you'll want to substract the areas below the zero line from the area above the zero line to integrate the area under the curve. When data below \code{thresh} make no sense for your question, you are able to set it. Then, all y-values below \code{thresh} are set to \code{thresh} and the regarding areas are not subtracted from the total area. }
  \item{dens }{ By default the scatter points are densified with factor 100. This means that the data density is increased by adding 100 data points between given adjacent data points by linear interpolation along x and y. When a threshold is set, this increases the accuracy of the result.}
  \item{sort.x }{ By default the vectors in \code{x} and \code{y} are ordered along increasing \code{x} because integration makes no sense with unordered data. You can override this by setting \code{sort.x = FALSE}
}
  \item{method }{
Specify how interpolation error should be estimated. Available methods include \code{"leave out"}, \code{"bootstrap"}, \code{"sorted bootstrap"}, \code{"constrained bootstrap"}, \code{"jackknife"}, \code{"jack-validate"}. True bootstrap is only effective when \code{sort.x = FALSE}. See details.
}
  \item{lo }{
When estimating interpolation error with \code{"leave out"} or \code{"jack-validate"}, how many data points should be left out randomly? Defaults to 2. See \code{method} and details.
}
  \item{it }{
How many iterations should be run when using \code{auc.mc} to estimate the interpolation error. Defaults to 100.
}
  \item{\dots }{
Any arguments passed through to \code{auc}.
}
}

\details{
When integrating the underlying assumption is that values can be interpolated linearly between adjacent data points. In many cases this is questionable. For estimating the linear interpolation error from the data at hand one may use Monte Carlo resampling methods. In \code{auc.mc} the following approaches are available:
\itemize{
  \item \code{leave out}: In each run \code{lo} data points are randomly omitted. This is quite straightforward, but the number of data points left out (\code{lo}) is arbitrary and thus the error terms estimated with this approach may be hardly defensible.
  \item \code{bootstrap}: Data are bootstrapped (sampling with replacement). Thus, some data points may repeat whereas others may be omitted. Due to the random sampling the order of data points is changed which may be unwanted with times series and may produce largely exaggerated error terms. This is only effective if \code{sort.x = FALSE}.
  \item \code{sorted bootstrap}: Same as before but ordering along \code{x} after bootstrapping may cure some problems of changed order. However, due to repeated data points time series spreading seasons but having data showing distinct seasonality may still be misrepresented.
  \item \code{constrained bootstrap}: Same as before but after ordering repeated data points are omitted. Thus, this equals leaving some measurements out at each run with a random number of leave outs. Numbers of leave outs typically show normal distribution around 3/4n.
  \item \code{jackknife}: \code{auc} is calculated for all possible combinations of \code{length(x)-1} data points. Depending on \code{length(x)} the number of combinations can be quite low.
  \item \code{jack-validate}: \code{auc} is calculated for all possible combinations of \code{length(x)-lo} to \code{length(x)-1} data points. Partly cures the "arbitrarity" problem of the \code{cross-validate} approach and produces stable summary statistics.
}
}

\value{
\code{auc} returns a numeric value that expresses the area under the curve. The unit depends from the input.

\code{auc.mc} returns a numeric vector containing the \code{auc} values of the \code{it} permutations. Just calculate summary statistics from this as you like. Due to the sampling approaches means and medians are not stable for most of the methods. \code{jackknife} and  \code{jack-validate} produce repeatable results, in the case of \code{cross-validate} it depends on n (\code{length(x)}) and \code{it}.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de}
}

\seealso{ \code{\link[caTools]{trapz}}, \code{\link[stats]{integrate}}}

\examples{
## Construct a data set (Imagine 2-hourly ghg emission data
## (methane) measured during a day).
## The emission vector (data in mg CH4 / m2*h) as a time series.
ghg <- ts(c(12.3, 14.7, 17.3, 13.2, 8.5, 7.7, 6.4, 3.2, 19.8, 
22.3, 24.7, 15.6, 17.4), start=0, end=24, frequency=0.5)
## Have a look at the emission development.
plot(ghg)
## Calculate what has been emitted that day
## Assuming that emissions develop linearly between
## measurements
auc(time(ghg), ghg)

## Test some of the auc.mc approaches
## "leave out" as default
auc.rep <- auc.mc(time(ghg), ghg)
## mean and median are well below the original value
summary(auc.rep)
## results for "bootstrap" are unstable (run several times)
auc.rep <- auc.mc(time(ghg), ghg, "boot")
summary(auc.rep)
## results for "jack-validate" are stable (run several times)
auc.rep <- auc.mc(time(ghg), ghg, "jack-val", lo=3)
summary(auc.rep)

## The effect of below.zero:
## Shift data, so that we have negative emissions (immissions)
ghg <- ghg-10
## See the difference
plot(ghg)
abline(h=0)
## With thresh = NULL the negative emissions are subtracted
## from the positive emissions
auc(time(ghg), ghg)
## With thresh = 0 the negative emissions are set to 0
## and only the emissions >= 0 are counted.
auc(time(ghg), ghg, thresh = 0)

}

\keyword{ math }
\keyword{ models }
\keyword{ ts }


