\name{export}
\alias{export}
\title{
simple export wrapper
}
\description{
Export your flux estimations easily
}

\usage{
export(x, digits = 4, ...)
}

\arguments{
  \item{x}{
A fluxes object.
}
  \item{digits}{
The number of digits that all numeric values in the output table shall have. Defaults to 4.
}
  \item{\dots}{
Further arguments to \code{\link[utils]{write.table}}. For instance to specify another field separator (defaults to tab delimited output files) and a file where to write the results to.
}
}

\details{
It's really very simple.
}

\value{
The function is called for its side effects. Nothing is returned.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\seealso{
\code{\link{flux}}, \code{\link{conz.prep}}, 
(also for examples)
}

\keyword{ manip }
