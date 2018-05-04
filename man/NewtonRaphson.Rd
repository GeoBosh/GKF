\name{NRUpdatingStep}
\alias{NRUpdatingStep}
\title{
Newton-Raphson updating step for computing the posterior mode.
}
\description{
  
  This function allows the computation of the posterior mode for
  nonlinear non-Gaussian observation model as described in
  \insertCite{@@jungbacker2007monte, chapter 3;textual}{GKF}. More
  precisely, we implement Eq(14) of
  \insertCite{jungbacker2007monte;textual}{GKF} following the details in
  (A10) and (A11).

}
\usage{
NRUpdatingStep(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,checkInputs = TRUE)
}
\arguments{
  \item{a0}{A \code{vector} giving the initial value/estimation of the
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
    observations. \dQuote{NA}-values are allowed (see \bold{Details}).}
  \item{checkInputs}{A \code{logical} stating whether the input shall be
    checked for consistency (\dQuote{storage.mode}, \dQuote{class}, and
    dimensionality, see \bold{Details}).}
}
\details{

  The function calls first \code{FKF} to compute the forward Kalman
  filter steps as described in \insertCite{@@jungbacker2007monte, A1}{GKF} 
  and then computes the backward recursions (A10) and (A11). As in
  \code{FKF}, \dQuote{NA}-values in the observation matrix \code{yt}
  (\code{xt} in \insertCite{jungbacker2007monte;textual}{GKF} notations)
  are supported.

}
\value{

  If \eqn{d}{d} is the dimension of the observations, and \eqn{n} the
  number of time observations, then the result is \eqn{d \times n}{d *
  n}-matrix containing the value of the updated posterior signal mode at
  each time step (in each column).

}
\references{
  \insertAllCited{}
}
\author{
  Tarak Kharrat
}
\note{

  Note the following notation difference with
  \insertCite{jungbacker2007monte;textual}{GKF}. Our \code{yt} is equivalent
  to \insertCite{jungbacker2007monte;textual}{GKF}'s \code{xt} and \code{Ht}
  is equivalent to \code{At}.

}

\seealso{
        \code{GaussianSignalSmoothing}
}
\examples{
"See Examples in package-vignette"
}
\keyword{Kalman-Filter}
\keyword{Posterior-Signal-mode}
