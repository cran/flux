\encoding{UTF-8}
\name{reco}
\alias{reco}
\title{
Model R_eco from CO2 exchange closed chamber data.
}
\description{
Model \eqn{R_{eco}}{R_eco} from CO2 exchange closed chamber data.
}

\usage{
reco(R, Temp, Tref = 10, T0 = -46.02, method = "all")
}

\arguments{
  \item{R}{
Numeric vector with ecosystem respiration (\eqn{R_{eco}}{R_eco}) flux rates.
}
  \item{Temp}{
Numeric vector with corresponding temperature values.
}
  \item{Tref}{
Numeric value giving the reference temperature used in the Arrhenius type model. Defaults to 10 (°C).
}
  \item{T0}{
Numeric value giving the activation temperature used in the Arrhenius type model. Defaults to -46.02 (°C).
}
  \item{method}{
Specifies the model to be used. Partial match is possible. One can either check all included models (\code{"all"}), select the best performing model (according to \code{\link[stats]{AIC}}; \code{"best"}), or do all models that do not fail in fitting (\code{"not.failed"}) or specify the wanted model directly: \code{"linear"}, \code{"arrhenius"}, \code{"Q10"}, \code{"lt"} (Lloyd & Taylor), \code{"ltr"} (Lloyd & Taylor restricted), \code{"logistic"}. See details.
}
}

\details{
Respiration is controlled by both biological and physical factors. Work by Arrhenius and van’t Hoff in the late-19th century on the temperature dependence of chemical reactions lead to the insight that there is a certain relationship between temperature and respiration (see review by Lloyd and Taylor, 1994). The most prominent models that have been used extenively in the literature can be fitted with this function. For an in-depth review, even more models and references see Richardson et al. 2006.

Models (T = \code{Temp}erature):
 \tabular{ll}{
    \code{linear}
    \tab \eqn{ R = \theta_1 + \theta_2*T }{ R = t1 + t2*T }
    \cr
    \code{arrhenius}
    \tab \eqn{ R = \theta_1 + exp\left[E_0 \left( \frac{1}{T_{Ref}-T_0} - \frac{1}{T-T_0} \right) \right]}{R = t1 * exp(E0 * (1/(Tref-T0) - 1/(T-T0)))}
    \cr
    \code{Q10}
    \tab \eqn{ R = \theta_1 \theta_2 ^{(T-T_{Ref})/10} }{ R = t1 * t2^((Temp-Tref)/10) }
    \cr
    \code{lt}
    \tab \eqn{ R = \theta_1 exp \left( \frac{-\theta_2}{T + 273.15 - \theta_3} \right) }{ R = t1 * exp(-t2 / (T+273.15-t3)) }
    \cr
    \code{ltr}
    \tab \eqn{ R = \theta_1 exp \left( \frac{-308.56}{T + 46.02} \right) }{ R = t1 * exp(-308.56 / (T + 46.02)) }
    \cr
    \code{logistic}
    \tab \eqn{ R = \frac{\theta_1}{1 + exp(\theta_2-\theta_3T)} }{ R = t1 / (1+exp(t2-t3*T)) }
    \cr    
 }
}

\value{
Either returns a list of \eqn{R_{eco}}{R_eco} models or the specified model structure. The wanted or resultant model can be fed into \code{\link{gpp}} or used on its own to predict Reco values.
}

\references{
Lloyd J, Taylor JA, 1994. On the temperature dependence of soil respiration. Functional Ecology 8:315–323.

Richardson et al. 2006. Comparing simple respiration models for eddy flux and dynamic chamber data. Agricultural and Forest Meteorology 141:219–234.
}

\author{
Gerald Jurasinski, \email{gerald.jurasinski@uni-rostock.de}
}

\note{
In its current implementation the \code{lt} and \code{logistic} models are easily over parameterized and therefore find singular gradients and provide no fit.
}

\seealso{
\code{\link{gpp}}
}

\examples{
## see examples at gpp
}

\keyword{ univar }
