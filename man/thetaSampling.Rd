\name{thetaSampling}
\alias{thetaSampling}
\alias{GaussianthetaSampling}

\title{
  Draw from the posterior distribution of the signal vectors
}
\description{
  The provided functions simulate M draw from the posterior distribution of the
  signal vector \eqn{P(\theta |Y_n)} in the linear Gaussian case as well
  as in the nonlinear non-Gaussian case.
}
\usage{
thetaSampling(a0, P0, dt, ct, Tt, Zt, Ht, Qt, yt, thetaHat, M,
              checkInputs = TRUE)
GaussianthetaSampling(a0, P0, dt, ct, Tt, Zt, Ht, Qt, yt, M,
              checkInputs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a0}{A \code{vector} giving the initial value/estimate of the state variable.}
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
    observations. \code{NA} values are allowed (see \bold{Details}).}
  \item{thetaHat}{A \code{matrix} containing the value of the signal
    posterior mode at each time step. Only required in the nonlinear
  non-Gaussian case (see \bold{Details}).}
  \item{M}{Number of samples required.}
  \item{checkInputs}{A \code{logical} stating whether the input shall be
    checked for consistency (\dQuote{storage.mode}, \dQuote{class}, and
    dimensionality, see \bold{Details}).}
}
\details{
  
  In the linear Gaussian case we follow the mean corrections algorithm
  described in \insertCite{@@durbin2012time, chapter 4;textual}{GKF}.
  The nonlinear non-Gaussian case follows
  \insertCite{@@jungbacker2007monte, Theorem 2;textual}{GKF}. Both
  functions make use of quantities computed by \code{FKF}.

}
\value{
  A \code{list} with the following elements:

  \item{\code{theta}}{A \eqn{d \times n \times M}{d*n*M}-array
    containing the \eqn{M} smoothed signal samples
    \eqn{\theta^i \sim f(\theta|Y_n), i=1,\dots,M}.
  }

  \item{logLik}{A \eqn{M}{M}-vector containing the log-likelihood
    \eqn{\log(P(\theta^i|Y_n))} associated to each sample.
  }
}

\references{
  \insertAllCited{}
}
\author{
  Tarak Kharrat and Georgi N. Boshnakov
}
\note{

  Note the following notation difference with
  \insertCite{@@jungbacker2007monte;textual}{GKF}: our \code{yt} is
  equivalent to \insertCite{@@jungbacker2007monte;textual}{GKF}'s
  \code{xt}, and \code{Ht} is equivalent to \code{At}.

}
\seealso{
  \code{\link{FKF}}, \code{\link{NRUpdatingStep}}
}
\examples{
## See Examples in the vignette:
vignette("vignette", package = "GKF")
}
\keyword{ts}
\keyword{Kalman}
