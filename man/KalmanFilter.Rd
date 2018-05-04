\name{FKF}
\alias{FKF}
\title{
  Fast Kalman Filter
}
\description{

  This function allows for fast and flexible Kalman filtering. Both, the
  measurement and transition equation may be multivariate and parameters
  are allowed to be time-varying. In addition, missing values,
  \code{NA}, in the observations are supported.

}
\usage{
FKF(a0, P0, dt, ct, Tt, Zt, Ht, Qt, yt, ComputeThetaLik = FALSE, 
    ThetaVal = list(), checkInputs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
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

  \item{ComputeThetaLik}{A \code{logical} stating whether the signal
    log-likelihood should be computed, see \bold{Details}). }

  \item{ThetaVal}{If \code{ComputeThetaLik} is \code{TRUE}, a
    \code{list} containing the signal sample where the log-likelihood
    should be evaluated, see \bold{Details}).}

  \item{checkInputs}{A \code{logical} stating whether the input shall be
    checked for consistency (\dQuote{storage.mode}, \dQuote{class}, and
    dimensionality, see \bold{Details}).}
}
\details{
  
  The following notation is closest to the one used in package
  \code{fkf} and similar to \insertCite{durbin2012time}{GKF}. The matrix
  \code{Ht} is not required to be definite positive anymore but just
  invertible. The computation is also very close to
  \insertCite{@@jungbacker2007monte, Theorem 1}{GKF}.

  \strong{State space form:}

  Let \eqn{m}{m} be the dimension of the state variable, \eqn{d}{d} be
  the dimension of the observations, and \eqn{n} the number of time
  observations. The transition equation and the measurement equation are given by
  \deqn{\alpha_{t + 1} = d_t + T_t \cdot \alpha_t + \eta_t}{
        alpha(t + 1) = d(t) + T(t) alpha(t) + eta(t)}
   \deqn{y_t = c_t + Z_t \cdot \alpha_t + \epsilon_t,}{
         y(t) = c(t) + Z(t) alpha(t) +  epsilon(t),}
  where \eqn{\eta_t}{eta(t)} and \eqn{\epsilon_t}{epsilon(t)} are iid
  \eqn{N(0, G_t)}{N(0, G(m))} and iid \eqn{N(0, H_t)}{N(0, H(d))},
  respectively, and \eqn{\alpha_t}{alpha(t)} denotes the state
  variable. We also define the signal \eqn{\theta_t}{theta(t)} by
  \deqn{\theta_t = c_t + Z_t \cdot \alpha_t} One can note that the
  signal has the dimension as the state.  The parameters admit the
  following dimensions:

  \tabular{lll}{
    \eqn{a_t \in R^m}{a[t] in R^m} \tab \eqn{d_t \in R^m}{d[t] in R^m} \tab \eqn{\eta_t \in R^m}{eta[t] in R^m} \cr
    \eqn{T_t \in R^{m \times m}}{d[t] in R^(m * m)} \tab \eqn{G_t \in R^{m \times m}}{d[t] in R^(m * m)} \tab \cr
    \eqn{y_t,\theta_t \in R^d}{y[t] in R^d} \tab \eqn{c_t \in R^d}{c[t] in R^d} \tab \eqn{\epsilon_t \in R^d}{epsilon[t] in R^d}. \cr
    \eqn{Z_t \in R^{d \times m}}{Z[t] in R^(d * m)} \tab \eqn{H_t \in R^{d \times d}}{H[t] in R^(d * d)} \tab
  }

\strong{NA-values:}

  \code{NA} values in the observation matrix \code{yt} are supported.
  If particular observations \code{yt[,i]} contain \code{NA}'s, the
  \code{NA} values are removed and the measurement equation is adjusted
  accordingly \insertCite{@see @durbin2012time, chapter 4.10}{GKF}.
  When the full vector \code{yt[,i]} is missing, the Kalman filter
  reduces to a prediction step.

% <------------------------------------->
\strong{Parameters:}

  The parameters can either be constant or deterministic
  time-varying. Assume the number of observations is \eqn{n}
  (i.e. \eqn{y = (y_t)_{t = 1, \ldots, n}, y_t = (y_{t1}, \ldots,
  y_{td})}{y = y[,1:n]}). Then, the parameters admit the following
  classes and dimensions:

  \tabular{ll}{
    \code{dt} \tab either a \eqn{m \times n}{m * n} (time-varying) or a \eqn{m \times 1}{m * 1} (constant) matrix. \cr
    \code{Tt} \tab either a \eqn{m \times m \times n}{m * m * n} or a \eqn{m \times m \times 1}{m * m * 1} array. \cr
    \code{Qt} \tab either a \eqn{m \times m \times n}{m * m * n} or a \eqn{m \times m \times 1}{m * m * 1} array. \cr
    \code{ct} \tab either a \eqn{d \times n}{d * n} or a \eqn{d \times 1}{d * 1} matrix. \cr
    \code{Zt} \tab either a \eqn{d \times m \times n}{d * m * n} or a \eqn{d \times m \times 1}{d * m * 1} array. \cr
    \code{Ht} \tab either a \eqn{d \times d \times n}{d * d * n} or a \eqn{d \times d \times 1}{d * d * 1} array. \cr
    \code{yt} \tab a \eqn{d \times n}{d * n} matrix.
  }

  If \code{ComputeThetaLik} is \code{TRUE}, the signal log-likelihood
  will be computed at the different \eqn{\theta^i}{} values contained in
  the \code{array} nammed "theta" (in \code{ThetaVal} which is a
  \code{list}).

  If \code{checkInputs} is \code{TRUE} each argument will be checked for
  correctness of the dimensionality, storage mode, and class. When
  possible, the dimension will be adapted by replication (a time
  constant matrix will be replicated for each time step internally).

}

\value{
  An S3 \code{list} with the following elements:

  \item{rl}{ \code{at} A \eqn{m \times (n + 1)}{m * (n + 1)}-matrix
    containing the predicted state variables, i.e. \eqn{a_t = E(\alpha_t
      | y_{t - 1})}{at[,t] = E(alpha[t] | y[,t- 1])}.}

  \item{Pt}{A \eqn{m \times m \times (n + 1)}{m * m * (n + 1)}-array
    containing the variances of \code{at}, i.e. \eqn{P_t = var(\alpha_t
      | y_{t - 1})}{Pt[,,t] = var(alpha[t] | y[,t - 1])}.}

  \item{att}{A \eqn{m \times n}{m * n}-matrix containing the filtered
    state variables, i.e. \eqn{a_{t|t} = E(\alpha_t | y_{t})}{at[,t] =
      E(alpha[t] | y[,t])}.}

  \item{Ptt}{A \eqn{m \times m \times n }{m * m * n}-array containing
    the variances of \code{att}, i.e. \eqn{P_{t|t} = var(\alpha_t |
      y_{t})}{Ptt[,,t] = var(alpha[t] | y[,t])}.}

  \item{vt}{A \eqn{d \times n}{d * n}-matrix of the prediction errors
    given by \eqn{v_t = y_t - c_t - Z_t a_t}{vt[,t] = yt[,t] - ct[,t] -
      Zt[,,t] \%*\% at[,t]}.}

  \item{Ft}{A \eqn{d \times d \times n}{d * d * n}-array which contains
    the variances of \code{vt}, i.e. \eqn{F_t = var(v_t)}{Ft[,,t] =
      var(v[,t])}.}

  \item{Fti}{A \eqn{d \times d \times n}{d*d*n}-array which contains the
    inverse of \code{Ft}, i.e. \eqn{F_t^i = F_t^{-1}}{Fti[,,t] =
      inv(Ft[,,t])}.}

  \item{Kt}{ A \eqn{m \times d \times n}{m * d * n}-array containing the
    \dQuote{Kalman gain} (see Jungbacker-Koopman(2007,A1)).}

  \item{logLik}{The log-likelihood \eqn{\log (P(Y_n))}.}

  \item{logLik_details}{The log-likelihood at each time step (useful for
    debug-purposes).}

  \item{thetalogLik}{The log-likelihood \eqn{\log (P(\theta^i))} for each
    signal sample \eqn{\theta^i,i=1,\dots,M}.}

  \item{thetalogLik_details}{The log-likelihood at each time step for
    each signal sample (useful for debug-purposes).}

  The first element of both \code{at} and \code{Pt} is filled with the
  function arguments \code{a0} and \code{P0}.

}
\references{
        % bibentry: durbin2012time
Durbin J and Koopman SJ (2012).
\emph{Time series analysis by state space methods},  number 38.
Oxford University Press.
% end:bibentry:  durbin2012time

        % bibentry: jungbacker2007monte
Jungbacker B and Koopman SJ (2007).
``Monte Carlo estimation for nonlinear non-Gaussian state space models.''
\emph{Biometrika}, \bold{94}(4), pp. 827--839.
% end:bibentry:  jungbacker2007monte

        % bibentry: sanderson2010armadillo
Sanderson C and others (2010).
``Armadillo: An open source C++ linear algebra library for fast prototyping and computationally intensive experiments.''
\emph{Report Version}, \bold{2}.
% end:bibentry:  sanderson2010armadillo
}
\author{
 Tarak Kharrat and Georgi N. Boshnakov
}
\note{

  Note the following notation difference with
  \insertCite{jungbacker2007monte;textual}{GKF}. Our \code{yt} is equivalent
  to \insertCite{jungbacker2007monte;textual}{GKF}'s \code{xt} and \code{Ht}
  is equivalent to \code{At}.
  
  The matrix inverse is computed by the \code{\link{InvMat}} routine (see its
  documentation for more details).
}
\examples{
"See Examples in package-vignette"
}
\keyword{Kalman-Filter}
