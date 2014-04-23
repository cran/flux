\encoding{UTF-8}
\name{round.POSIXlt}
\alias{round.POSIXlt}
\alias{round.POSIXct}
\title{
Round times.
}
\description{
There are \code{\link[base]{round}} methods in base for objects of \code{\link[base]{DateTimeClasses}}. However, they can only round to full second, minutes, hours, and days. These functions offer some more options.
}

\usage{
\method{round}{POSIXlt}(x, digits = c("mins", "5mins", "10mins", "15mins", "quarter hours", 
"30mins", "half hours", "hours"))

\method{round}{POSIXct}(x, ...)
}

\arguments{
  \item{x}{
A \code{\link[base]{DateTimeClasses}} object to be rounded. Needs to carry a timezone specification.
}
  \item{digits}{
Either a character string specifying the time units to round to (see choices above) or a numeric specifying the minutes to round to. To go to seconds just use values < 1, to go beyond the hour just use values > 60.
}
  \item{\dots}{
Further arguments to methods.
}
}

\value{
A \code{POSIXct} object with rounded, not truncated date times.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de},

borrowing heavily from https://stat.ethz.ch/pipermail/r-help/2012-June/315336.html
}

\examples{
# Current time in GMT and as class "POSIXlt"
zlt <- as.POSIXlt(Sys.time(), "GMT") 

# Same time as class POSIXct
zct <- as.POSIXct(zlt)

# round to minute
round(zct)

# round to half hour
round(zct, "30mins")
round(zct, "half hour")
round(zct, 30)

# round to 20 minutes
round(zlt, 20)

# round to 30 seconds
round(zlt, 0.5)
}

\keyword{ manip }
