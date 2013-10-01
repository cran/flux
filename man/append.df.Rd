\name{append.df}
\alias{append.df}
\title{
Append a data.frame to another including consistency checks
}
\description{
Often ghg concentration data come in chunks. This function provides a wrapper for appending data.
}

\usage{
append.df(orig, add)
}

\arguments{
  \item{orig}{
A data.frame
}
  \item{add}{
Another data.frame
}
}

\details{
The two data.frames are appended based on common columns. A warning is issued if some column names do not match. New columns are silently added.
}

\value{
Data.frame 
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\examples{
## add later
}

\keyword{ manip }
