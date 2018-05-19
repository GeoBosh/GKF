\name{PosteriorSignalMode}
\alias{PosteriorSignalMode}

\title{
        Posterior Signal mode computation
}
\description{
        The function iterates the Newton-Raphson updating step implemented in 
        \code{NRUpdatingStep} until convergence to obtain the posterior signal mode.
}
\usage{
PosteriorSignalMode(a0, P0, dt, ct, Tt, Zt, Qt, yt, jac, hess, 
                    Pars, g = 0, tol = 1e-06, maxIter = 10, 
                    checkInputs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{a0}{A \code{vector} giving the initial value/estimation of the state variable.}
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
    observations. \dQuote{NA}-values are allowed (see \bold{Details}).}
  \item{jac}{ A function to compute the jacobian \eqn{\frac{\partial p(y|\theta)}{\partial \theta}} 
  (see \bold{Details}.)}
  \item{hess}{
  A function to compute the hessian \eqn{\frac{\partial^2 p(y|\theta)}{\partial \theta^2}}, 
  (see \bold{Details}).
}
  \item{Pars}{ A \code{list} containing distribution specific parameters. Used by \code{jac} and 
               \code{hess} (see \bold{Details}). 
}
  \item{g}{ A \code{scalar}, \code{vector} or a \code{matrix} containing the signal initial guess.
   (see \bold{Details}).
}
  \item{tol}{ Convegence is declared if the mean squared error between two succesive iteration is smaller
  than \code{tol}, default set to 1e-6 .  
}
  \item{maxIter}{ the iterations stop if the maximum number of iteration \code{maxIter} is reached, default 10.
}
  \item{checkInputs}{A \code{logical} stating whether the input shall be
    checked for consistency (\dQuote{storage.mode}, \dQuote{class}, and
    dimensionality, see \bold{Details}).}
}
\details{
  
  The function implements the Newton-Raphson iterations suggested by
  \insertCite{jungbacker2007monte;textual}{GKF}.
    
  \strong{Jacobian and Hessian computation:}
     
  In order to compute Eq. (14) from
  \insertCite{jungbacker2007monte;textual}{GKF} one needs to obtain the
  jacobian and the hessian of \eqn{\log p(y_t|\theta_t)} at different
  points \eqn{y_t} and \eqn{\theta_t}. Hence, user needs to provide
  \code{R} routines to do so. We require those function to have the
  following signatures:
  \tabular{rl}{
    \code{jac} \tab \code{function(theta,y,Pars)} and returns a \eqn{d \times 1} vector. \cr

    \code{hess} \tab  \code{function(theta,y,Pars)} and returns a \eqn{d \times d} matrix. \cr 
  }
  \code{Pars} is a \code{list} with some extra parameters used by
  \code{jac} and \code{hess}. If no extra parameters are needed, an
  empty \code{list} should be provided.  The functions \code{jac} and
  \code{hess} should handle \code{NA} (missing values) in \eqn{yt} and
  still return \code{numeric} vectors or matrices; \code{NA} values in
  \eqn{\theta_t} are not expected.

}
\value{
  A \code{list} with the following elements:

  \item{\code{theta_hat}}{If \eqn{d}{d} is the dimension of the
    observations, and \eqn{n} the number of time observations,
    \code{theta_hat} is \eqn{d \times n}{d *n}-matrix containing the
    value of the updated posterior signal mode at each time step (in
    each column).}

  \item{\code{Iter}}{Number of iteration to reach convergence.}
  \item{\code{err}}{value of the mean squared error at the final step.}
}
\references{
  \insertAllCited{}
}
\author{
  Tarak Kharrat and Georgi N. Boshnakov
}
\seealso{
  \code{\link{NRUpdatingStep}}, \code{\link{GaussianSignalSmoothing}}
}
\examples{
## See Examples in the vignette:
vignette("vignette", package = "GKF")
}
\keyword{ts}
\keyword{Kalman}
