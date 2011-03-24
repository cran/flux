\encoding{UTF-8}
\name{plot.gpp}
\alias{plot.gpp}
\title{
Plot diagnostic plots for GPP (NEE) models derived with reco and gpp.
}
\description{
Plot diagnostic plots for GPP (NEE) models derived with reco and gpp.
}

\usage{
\method{plot}{gpp}(x, nm = "", ...)
}

\arguments{
  \item{x}{
Object of class \code{gpp} that is returned by \code{\link{gpp}}.
}
  \item{nm}{
The three panels of the resulting plot are already be named. However, if you'd like to add something you can do it here.
}
  \item{\dots}{
Further arguments passed to \code{\link[graphics]{plot.default}}.
}
}

\details{
The function produces a three panel plot representing in this order from left to right: (1) \eqn{R_{eco}}{R_eco} plot and the used \eqn{R_{eco}}{R_eco} model. (2) Combined plot of the NEE/GPP data with the measured NEE vs PAR, the derived GPP and the modelled \eqn{R_{eco}}{R_eco}. (3) Diagnostic plot of \eqn{NEE_{measured}}{R_measured} vs \eqn{NEE_{predicted}}{R_predicted}
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

\keyword{ hplot }