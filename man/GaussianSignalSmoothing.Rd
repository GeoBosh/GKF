\name{GaussianSignalSmoothing}
\alias{GaussianSignalSmoothing}
\title{
  Gaussian signal smoothing.
}
\description{
  The function apply Kalman smoother to compute smoothed values of the
  signal vectors, together with their variance-covariance matrices in
  the linear Gaussian case.
}
\usage{
GaussianSignalSmoothing(a0, P0, dt, ct, Tt, Zt, Ht, Qt, yt,  
                        checkInputs = TRUE)
}
\arguments{
  \item{a0}{A \code{vector} giving the initial value/estimate of the
    state variable.} 
  \item{P0}{A \code{matrix} giving the variance of \code{a0}.}
  \item{dt}{A \code{matrix} giving the intercept of the transition
    equation (see \bold{Details}).}
  \item{ct}{A \code{matrix} giving the intercept of the measurement
    equation (see \bold{Details}).}
  \item{Tt}{An \code{array} giving the factor of the transition
    equation (see \bold{Details}).}
  \item{Zt}{An \code{array} giving the factor of the measurement
    equation (see \bold{Details}).}
  \item{Qt}{An \code{array} giving the variance of the innovations of
    the transition equation (see \bold{Details}).}
  \item{Ht}{An \code{array} giving the variance of the disturbances of
    the measurement equation (see \bold{Details}).}
  \item{yt}{A \code{matrix} containing the
    observations. \\code{NA}-values are allowed (see \bold{Details}).}
  \item{checkInputs}{A \code{logical} stating whether the input shall be
    checked for consistency (\dQuote{storage.mode}, \dQuote{class}, and
    dimensionality, see \bold{Details}).}
}
\details{

  The computation follows the Durbin-Koopman(2012,chap 4.5.3)
  \insertCite{@@durbin2012time, chapter 4.5.3}{GKF}, backward recursions
  and calls \code{FKF} to compute the Kalman Filter quantites
  beforehand.

}
\value{
  A \code{list} with the following elements:

  \item{theta_hat}{A \eqn{d \times n}{d*n} matrix containing the
    smoothed signal, i.e.  \eqn{\hat{\theta}= E(\theta|Y_n)}.}

  \item{Vtheta_hat}{A \eqn{d \times d \times n}{d*d*n} array containing
    the variance/covariance matrices , i.e
    \eqn{V_t=var(\theta_t|Y_n),t=1,\dots,n}.}
}
\references{
  \insertAllCited{}
}
\author{
  Tarak Kharrat and Georgi N. Boshnakov
}
\seealso{
\code{NRUpdatingStep}
}
\examples{
## See Examples in the vignette
}
\keyword{Signal-Smoothing}
\keyword{Posterior-Signal-mode}
