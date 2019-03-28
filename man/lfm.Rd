\name{lfm}
\alias{lfm}
\title{Estimate the linear feedback model in Blundell, Griffith and Windmeijer
  (2002)}
\description{Estimate the linear feedback model in Blundell, Griffith and
  Windmeijer "Individual effects and dynamics in count data models", Journal of
  Econometrics 108 (2002) 113-131}
\usage{
lfm(formula, data, effect = "individual", model = "onestep",
    weight.matrix = "instruments", index = NULL, start = NULL)
}
\arguments{
\item{formula}{Similar to the \code{pgmm()} function in package \code{plm}.  A
  symbolic description for the model to be estimated.  Indicate a multi-part
  formula, the first two parts describing the covariates and the gmm
  instruments and, if any, the third part the 'normal' instruments. The first
  independent variable must be the lag of the dependent variable.}
\item{data}{A \code{pdata.frame}, or a \code{data.frame} if using the \code{index} option.}
\item{effect}{Either \code{"individual"} or \code{"twoways"}. The former only
  includes individual fixed effects while the latter also includes time fixed
  effects.}
\item{effect}{Either \code{"onestep"} or \code{"twosteps"}. Whether to do
  one-step GMM or two-step GMM.}
\item{weight.matrix}{Either \code{"identity"} or \code{"instruments"}. Whether
  to use the identity matrix of the cross product of the instruments for the
  first-step weight matrix.}
\item{index}{If \code{data} is not a \code{pdata.frame}, use this option to
        declare the names of the individual and time indexes.}
\item{start}{An optional vector of starting values for optimization.}
}
\value{
  \item{call}{The matched call}
  \item{coefficients}{The estimated coefficient}
  \item{fitted.values}{\code{data.frame} of fitted values}
  \item{first}{The first stage estimates}
  \item{fixed.effects}{Estimates of the individual fixed effects}
  \item{model}{The variables used for estimation for each individual}
  \item{residuals}{\code{data.frame} of residuals}
  \item{vcov}{The covariance matrix of the coefficients}
  \item{W1}{The first-stage weight matrix used}
  \item{W2}{The second-stage (efficient) weight matrix used (only returned if
    \code{model = "twosteps"} is used)}
  \item{Z}{The instrument matrix for each individual}
  \item{obj}{The value of the GMM objective function at the estimated parameters.}
}
\examples{
\dontrun{
lfm(y ~ lag(y, k = 1) + x | lag(y, k = 2:4) + lag(x, k = 1:4),
    data = data, effect = "individual", model = "onestep", index = c("i", "t"))
}
}
