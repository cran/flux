\encoding{UTF-8}
\name{lips}
\alias{lips}
\title{
Linear interpolation between data points similar to approx.
}
\description{
Linear interpolation between data points similar to approx. x may be a time vector.
}

\usage{
lips(x, y, x.step = 1)
}

\arguments{
  \item{x}{
Numeric vector or one that could be coerced to numeric along which interpolation shall take place. May be a time vector (POSIXlt or POSIXct).
}
  \item{y}{
Numeric vector of values which shall be interpolated.
}
  \item{x.step}{
Numeric giving at what time interval interpolation shall be done. In seconds! Thus half hourly interpolation is achieved with \code{x.step = 1800}.
}
}

\value{
Data.frame containing the interpolated x (\code{x.out}) and y (\code{x.out}) values.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de}
}

\seealso{
\code{\link[stats]{approx}}
}

\examples{
## has to be added
}

\keyword{ manip }
\keyword{ ts }
