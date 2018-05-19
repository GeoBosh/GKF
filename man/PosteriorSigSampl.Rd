\name{PosteriorSignalSampling}
\alias{PosteriorSignalSampling}

\title{
     Wrapper for Signal Sampling Smoother routines
}
\description{

  The function provides a common interface to simulate from the
  posterior distribution of the signal vector \eqn{P(\theta |Y_n)} in
  the linear Gaussian case, as well as in the nonlinear non-Gaussian
  case by wrapping \code{thetaSampling} and
  \code{GaussianthetaSampling}.

}
\usage{
PosteriorSignalSampling(a0, P0, dt, ct, Tt, Zt, Qt, yt, SpecificPars, M = 1, 
                        type = c("LinGaussian", "NLnonGaussian"))
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
  \item{yt}{A \code{matrix} containing the
    observations. \code{NA}-values are allowed (see \bold{Details}).}
  \item{SpecificPars}{A \code{list} containing the model specific parameters
    (see \bold{Details}).}
  \item{M}{Number of samples required.}   
  \item{type}{ Choose between \code{LinGaussian} for linear Gaussian 
    and \code{NLnonGaussian} for nonlinear non-Gaussian models.}
}
\details{

  \code{SpecificPars} contains the model specific parameters which are
  as follows:
        
  \strong{Linear Gaussian model:}
        
  The only required parameter is the \code{array} (\code{matrix})
  \code{Ht} giving the variance of the disturbances of the measurement
  equation.
        
  \strong{Nonlinear non-Gaussian model:}
        
  Routines to compute the jacobian (\code{jac}) and the hessian
  (\code{hess}) with respect to the signal are required. They should
  handle missing values (\code{NA}) in the observation vector \code{yt}
  and always return a \code{numeric} vector or a matrix. Otherwise,
  related computation will break. Other parameters can also be provided
  (see \code{PosteriorSignalMode}). If not, they will be set to their
  default values.

}
\value{
  A \code{list} with the following elements:
 
  \item{\code{theta}}{An \eqn{d \times n \times M}{d*n*M}-array
    containing the \eqn{M} smoothed signal samples
    \eqn{\theta^i \sim f(\theta|Y_n), i=1,\dots,M}.}

  \item{\code{CondthetaLogLik}}{A \eqn{M}{M}-vector containing the
    log-likelihood \eqn{\log (P(\theta^i|Y_n))} associated to each
    sample.}

  \item{\code{UncondobsLogLik}}{A \code{double} value representing the
    unconditional log-likelihood of the observation vector \eqn{\log
      (P(Y_n))}.}

  \item{\code{UncondthetaLogLik}}{An \eqn{M}{M}-vector containing the
    unconditional log-likelihood \eqn{\log (P(\theta^i))} associated
    with each sample.}

}
\references{
        % bibentry: jungbacker2007monte
Jungbacker B and Koopman SJ (2007).
``Monte Carlo estimation for nonlinear non-Gaussian state space models.''
\emph{Biometrika}, \bold{94}(4), pp. 827--839.
% end:bibentry:  jungbacker2007monte

        % bibentry: durbin2012time
Durbin J and Koopman SJ (2012).
\emph{Time series analysis by state space methods},  number 38.
Oxford University Press.
% end:bibentry:  durbin2012time

}
\author{
  Tarak Kharrat
}
\seealso{
  \code{\link{thetaSampling}},
  \code{\link{GaussianthetaSampling}},
  \code{\link{PosteriorSignalMode}}
}
\examples{      
## See Examples in the vignette:
vignette("vignette", package = "GKF")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{Kalman}

