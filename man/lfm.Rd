\name{lfm}
\alias{lfm}
\title{Estimate the linear feedback model in Blundell, Griffith and Windmeijer (2002)}
\description{Estimate the linear feedback model in Blundell, Griffith and Windmeijer "Individual effects and dynamics in count data models", Journal of Econometrics 108 (2002) 113-131}
\usage{
lfm(formula, data, effect = "individual", model = "onestep")
}
\arguments{
\item{formula}{Similar to the \code{pgmm()} function in package \code{plm}.  A
  symbolic description for the model to be estimated.  Indicate a multi-part
  formula, the first two parts describing the covariates and the gmm
  instruments and, if any, the third part the 'normal' instruments. The first
  independent variable must be the lag of the dependent variable.}
\item{data}{A \code{pdata.frame}.}
\item{effect}{Either \code{"individual"} or \code{"twoways"}. The former only includes individual fixed effects while the latter also includes time fixed effects.}
\item{effect}{Either \code{"onestep"} or \code{"twosteps"}. Whether to do one-step GMM or two-step GMM.}
}
\value{
  \item{call}{The matched call}
  \item{coefficients}{The estimated coefficient}
  \item{D}{The average of the Jacobian of the sample moment conditions over
    each individual}
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
}
\examples{
# Not run:
lfm(y ~ lag(y, k = 1) + x | lag(y, k = 2:3) + lag(x, k = 2:3), data = df)
}

