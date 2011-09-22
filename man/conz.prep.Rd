\encoding{UTF-8}
\name{conz.prep}
\alias{conz.prep}
\title{
Prepare data for flux rate estimation with \code{\link{flux}}
}
\description{
The function simply constructs a list of \code{\link{data.frame}}s that each contains the data for one dynamic closed chamber measurement. 
}

\usage{
conz.prep(dat, columns, factors, nmes, min.cm = 3)
}

\arguments{
  \item{dat}{
\code{\link{data.frame}} containing static closed chamber data for several chamber measurements. See columns and example for details.
}
  \item{columns}{
A character vector giving the names of the columns in \code{dat} that shall be considered. Defaults to NULL which means that all columns are taken.
}
  \item{factors}{
A character vector giving the names of the columns that are used to partition the data in \code{dat} into chunks that each contain the data for one chamber placement. \code{factors} \eqn{\subset} \code{columns}! See example.
}
  \item{nmes}{
A character vector giving the names of the columns that are used to name the data chunks. \code{nmes} \eqn{\subset} \code{columns}! See example.
}
  \item{min.cm}{
Integer giving the minimum number of concentration measurements allowed in the tables that are created per chamber measurement. Defaults to 3 because   a linear fit to 2 points does not make any sense. Attention: Chamber placements with less than \code{min.cm} measurements are quietly skipped. 
}
}

\details{
This could easily be hand scripted but the function shall provide a simple way to obtain the structure needed for \code{\link{flux}} and it also carries naming information.
}

\value{
Returns a list with 2 entries. The first is itself a list of \code{\link{data.frame}}s containing the concentration measurements that result from the field sampling during one chamber placement (if \code{factors} has been specified correctly) and the columns specified in \code{columns}. The entries in the list are named according to \code{nmes}. However, the second part of the upper level list is a table with the naming information. This is handed over to \code{\link{flux}} and \code{\link{plot.fluss}}. See example.
}

\author{
Gerald Jurasinski <gerald.jurasinski@uni-rostock.de>
}

\seealso{
\code{\link{flux}}
}

\examples{
## load example data
data(tt.pre)

## extract field concentration measurements
gcd <- tt.pre[tt.pre$sampletype_a=="P",]

## partition the data into data tables per chamber measurement
# first tell which columns of the data you want to carry on
columns <- c("date", "veg", "spot", "time_min", "temp_dC", "cham_area", 
"cham_vol", "CO2ppm", "CO2Code", "CH4ppb", "CH4Code", "N2Oppb", "N2OCode")
# then do the partitioning
gcd.parts <- conz.prep(gcd, columns, factors = c("date", "spot", "veg"), 
nmes = c("date", "veg", "spot"))
# have a look at the first three tables
gcd.parts$tables[1:3]
# have a look at the names part of the returned object
gcd.parts$nmes
# use inspect to have a look at (a) specific data table(s)
inspect(gcd.parts, c("2011-03-15.c.3", "2011-03-15.c.6", "2011-03-15.p.6"))
# inspect the same tables using their indices
inspect(gcd.parts, c(3,6,12))
inspect(gcd.parts, c("c.3", "c.6", "p.6"))
}

\keyword{ manip }
