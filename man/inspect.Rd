\name{inspect}
\alias{inspect}
\title{
Inspect and alter prepared ghg concentration data
}
\description{
The function allows straightforward inspection and alteration of ghg concentration data that have been prepared using \code{\link{conz.prep}}
}
\usage{
inspect(x, what, sustain = NULL, retain = FALSE)
}

\arguments{
  \item{x}{
Object returned by \code{\link{conz.prep}} containing tables with chamber measurement data.
}
  \item{what}{
Specifies the concentration measurement tables in \code{x} that ought to be inspected/altered. For inspection either give a numeric vector with the indices of the tables or a character vector with the names of the entries in \code{x}. For altering data tables (i.e. deleting some values) use a named list with the names the names of the tables and the entries either numeric values or vectors that specify the concentration measurements that shall be skipped.
}
  \item{sustain}{
The names of the data tables in \code{x} can be quite complex. However, a certain part may be the same for all the tables you want to inspect/alter. In this case just give this part here and skip it in \code{what}.
}
  \item{retain}{
Logical. When you alter \code{x} by giving a list structure to \code{what}, shall the original data tables be returned as well? Defaults to FALSE.
}
}

\details{
The typical workflow is to look at the diagnostic plots of the results of the \code{\link{flux}} estimation and then turn to \code{inspect} for having a closer look at the data or to delete some concentration measurements for further estimations.
}

\value{
Either the data tables to inspect are returned in a list or the altered \code{x}. In case \code{retain = TRUE} the original tables are appended to \code{x}. The respective list item is \code{tables.orig}.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\seealso{
\code{\link{flux}}, \code{\link{fluxes}}, \code{\link{conz.prep}}, 
(also for examples)
}

\keyword{ manip }
