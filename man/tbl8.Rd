\encoding{UTF-8}
\name{tbl8}
\alias{tbl8}
\title{
Extract the relevant data of a bulk model object to a data.frame
}
\description{
Extract the relevant data of an object of class "\code{breco}" to a \code{data.frame}
} 

\usage{
tbl8(models)
}

\arguments{
  \item{models}{
An object of class "\code{breco}" returned by \code{reco.bulk} or an object of class "\code{bgpp}" returned by \code{gpp.bulk}. 
}
}

\value{
Returns a \code{data.frame} with the model coefficients, brute force R2s (not reliable because non linear responses are fitted), timestamps and the \code{which.Temp} string or the \code{offset} for \eqn{R_{eco}}{R_eco} models and GPP models, respectively.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de},
based on an idea by Sascha Beetz, \email{sascha.beetz@uni-rostock.de}
}

\seealso{
\code{\link{reco}}, \code{\link{reco.bulk}}, \code{\link{gpp}}, \code{\link{gpp.bulk}}
}

\examples{
## see examples at reco.bulk
}

\keyword{ manip }
