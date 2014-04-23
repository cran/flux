\encoding{UTF-8}
\name{checkm}
\alias{checkm}
\title{
Bring modelled and measured values together based on timestamp
}
\description{
Trivial function that is a simple wrapper for frequent task: Bringing together the measured and modelled values, for instance to do a posteriori analyses of model performance. 
}

\usage{
checkm(modelled, measured, t.unit = NULL)
}

\arguments{
  \item{modelled}{
A \code{data.frame} with the modelled fluxes and all additional data reported. For instance resulting from doing \code{\link{budget.reco}} on a "\code{breco}" object. One variable has to be the timestamp and it should be named exactly like this.
}
  \item{measured}{
A \code{data.frame} with the originally measured fluxes and additional data. One variable has to be the timestamp and it should be named exactly like this.
}
  \item{t.unit}{
If NULL, data in \code{modelled} and \code{measured} are merged based on their timestamp by calculating the minimum difference in time between all entries (see Details). If !NULL, character string specifying to which time interval the timestamps in \code{measured} should be rounded using \code{\link{round.POSIXlt}}. May be one of \code{"mins"}, \code{"5mins"}, \code{"10mins"}, \code{"15mins"}, \code{"quarter hours"}, \code{"30mins"}, \code{"half hours"}, \code{"hours"}. Alternatively, a numeric specifying the minutes to round to. To go to seconds just use values < 1, to go beyond the hour just use values > 60.
}
}

\details{
\bold{Case 1 (t.unit = NULL)} Data are merged by calculating the difference in time between all timestamps in \code{modelled} and all timestamps in \code{measured} and identifying the minimum difference to each measured flux. If minimum difference between measured and modelled flux > 1h, no modelled flux is assigned. This approach is a bit slower but it is not necessary to give a correct \code{t.unit}, which makes it less error prone.

\bold{Case 2 (t.unit != NULL)} After rounding the timestamps in \code{measured} according to \code{t.unit} and transforming both timestamnps to character vectors \code{modelled} and \code{measured} are merged based on these  timestamps and only data rows that are present in both are retained. Therefore \code{t.unit} has to be specified according to the interval of the timestamps in \code{modelled}.
}

\value{
Data.frame containing the corresponding rows of \code{modelled} and \code{measured}
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de}
}

\examples{
## See examples at reco.bulk
}

\keyword{ manip }