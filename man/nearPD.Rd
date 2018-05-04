\name{nearPD}
\alias{nearPD}
\title{Nearest Positive Definite Matrix}
\description{
  Compute the nearest positive definite matrix to an approximate one,
  typically a correlation or variance-covariance matrix.
}
\usage{
nearPD(x,w, corr = FALSE, keepDiag = FALSE, EnforcePosDef = TRUE,
       doSym = FALSE, ensureSymmetry = FALSE,
       eig_tol = 1e-6, conv_tol = 1e-7, posd_tol = 1e-3,
       maxit = 100, chol = TRUE)
}
\arguments{
  \item{x}{numeric \eqn{n \times n}{n * n} approximately positive
    definite matrix, typically an approximation to a correlation or
    covariance matrix.  If \code{x} is not symmetric (and
    \code{ensureSymmetry} is not false), \code{\link{symmpart}(x)} is used.}
  \item{w}{ vector defining a diagonal weight matrix \code{diag(w)}. 
  If missing, a vector of ones is used.}
  \item{corr}{logical indicating if the matrix should be a
    \emph{correlation} matrix.}
  \item{keepDiag}{logical, generalizing \code{corr}: if \code{TRUE}, the
    resulting matrix should have the same diagonal
    (\code{\link{diag}(x)}) as the input matrix.}
  \item{EnforcePosDef}{logical indicating if a
    \code{\link[sfsmisc]{posdefify}()} eigen step should be applied to
    the result of the Higham algorithm.}
  \item{doSym}{logical indicating if \code{X <- (X + t(X))/2} should be
    done, \cr 
    after \code{X <- tcrossprod(Qd, Q)}; some doubt if this is necessary.}
  \item{ensureSymmetry}{logical; if true, \code{doSym} will be called.}
  \item{eig_tol}{defines relative positiveness of eigenvalues compared
    to largest one, \eqn{\lambda_1}. Eigen values \eqn{\lambda_k} are
    treated as if zero when \eqn{\lambda_k / \lambda_1 \le eig_tol}.}
  \item{conv_tol}{convergence tolerance for Higham algorithm.}
  \item{posd_tol}{tolerance for enforcing positive definiteness (in the
    final \code{posdefify} step when \code{do2eigen} is \code{TRUE}).}
  \item{maxit}{maximum number of iterations allowed.}
  \item{chol}{if true, the square root (Cholesky) decomposition of the
    resulting matrix.}   
  }
  \details{
    
  This implements the algorithm of
  \insertCite{higham2002computing;textual}{GKF} and then (if
  \code{EnforcePosDef} is true) forces positive definiteness using code
  from \code{\link[sfsmisc]{posdefify}}.  The algorithm of
  \insertCite{knol1989least;textual}{GKF} (not implemented here) is more
  general in (1) that it allows constraints to fix some rows (and
  columns) of the matrix and (2) to force the smallest eigenvalue to
  have a certain value.

  Note that setting \code{corr = TRUE} just sets \code{diag(.) <- 1}
  within the algorithm.

  \insertCite{higham2002computing;textual}{GKF}
  uses Dykstra's correction which is also the choice made by \code{nearPD}.
  
  if \code{chol} is true, the approximate Cholesky decomposition using
  filtered eigenvalues is computed and returned
  \insertCite{@see also @cheng1998modified}{GKF}.

}
\value{
  a list with components
  \item{mat}{a matrix: the computed positive-definite matrix.}
  \item{eigenvalues}{numeric vector of eigen values of \code{mat}.}
  \item{corr}{logical, just the argument \code{corr}.}
  \item{normF}{the Frobenius norm (\code{\link{norm}(x-X, "F")}) of the
  difference between the original and the resulting matrix.}
  \item{iterations}{number of iterations needed.}
  \item{converged}{logical indicating if iterations converged.}
  \item{CholDec}{the square root (Cholesky) decomposition of \code{mat}.}
}

\references{%% more in /u/maechler/R/Pkgs/sfsmisc/man/posdefify.Rd
  \insertAllCited{}
}
\author{Tarak Kharrat inspired by Higham matlab code}

\examples{
## Higham(2002), p.334f - simple example
A <- matrix(1, 3,3); A[1,3] <- A[3,1] <- 0
n.A_Rcpp <- nearPD(A, corr=TRUE)
stopifnot(all.equal(0.7593000,n.A_Rcpp$mat[1,2],tolerance=1e-6),
          all.equal(0.1577192,n.A_Rcpp$mat[1,3],tolerance=1e-6))
    
set.seed(27)
m <- matrix(round(rnorm(25),2), 5, 5)
m <- m + t(m)
diag(m) <- pmax(0, diag(m)) + 1
(m <- round(cov2cor(m), 2))
    
near.m_Rcpp <- nearPD(m)
stopifnot(all.equal(near.m_Rcpp$normF,0.9412029,tolerance=1e-6))
    
## A longer example, extended from Jens' original,
## showing the effects of some of the options:
pr <- matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
               0.477, 1,     0.516, 0.233, 0.682, 0.75,
               0.644, 0.516, 1,     0.599, 0.581, 0.742,
               0.478, 0.233, 0.599, 1,     0.741, 0.8,
               0.651, 0.682, 0.581, 0.741, 1,     0.798,
               0.826, 0.75,  0.742, 0.8,   0.798, 1),
             nrow = 6, ncol = 6)
    
nc1_Rcpp <- nearPD(pr, conv_tol = 1e-7,posd_tol  = 1e-4)
stopifnot(all.equal(nc1_Rcpp$normF,0.06304718,tolerance=1e-6),
          all.equal(min(nc1_Rcpp$eigenvalues),0.0004212751,tolerance=1e-6))
}
\keyword{Matrix-Factorization}
