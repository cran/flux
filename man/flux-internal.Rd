\encoding{UTF-8}
\name{flux-internal}
\alias{print.fluss}
\alias{print.fluxes}
\alias{print.fluxxes}
\title{ Internal functions of the package flux that are usually not called (or only implicitly) by the user }
\description{
Mainly print functions for giving nice output.
}

\usage{
\method{print}{fluss}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{fluxes}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{fluxxes}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{ object to be printed }
  \item{digits}{ number of digits }
  \item{\dots}{ Further arguments. }
}

\details{
    These are just printing methods to give nice output. Usually they are not to be called by the user.
}

\value{
  prints the output of a function and invisibly gives back the original object
}

\author{ Gerald Jurasinski }

\keyword{ internal }
