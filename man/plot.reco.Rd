\encoding{UTF-8}
\name{plot.reco}
\alias{plot.reco}
\title{
Plot diagnostic plots for Reco models.
}
\description{
Plot diagnostic plots for \eqn{R_{eco}}{R_eco} models.
}

\usage{
\method{plot}{reco}(x, ...)
}

\arguments{
  \item{x}{
Object of class \code{reco} that is returned by \code{\link{reco}}.
}
  \item{\dots}{
Further arguments passed to \code{\link[graphics]{plot.default}}.
}
}

\details{
The function produces a plot of the reco model.
}

\value{
The function is invoked for its side effects and does not return anything.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\seealso{
\code{\link{chop}}, \code{\link{fluxx}}, \code{\link{gpp}}, \code{\link{reco}}
}

\examples{
## see examples at gpp
}

\keyword{ hplot }