\name{Robust_chol}
\alias{Robust_chol}
\title{
     Modified Cholesky factorization
}
\description{

  This function computes an altenative Cholesky factorization more
  adapted to symmetric semi-definte matrices.  It is based on
  eigenvalues decomposition (equivalent to the singular value
  decomposition in the symmetric case) and shrinkage.

}
\usage{
Robust_chol(B, eig_tol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{A symmetric positive semi-definite \code{matrix} to be factorized.}
  \item{eig_tol}{
    defines relative positiveness of eigenvalues compared to largest
    one, \eqn{\lambda_1}.  Eigenvalues \eqn{\lambda_k} are treated as if
    zero when \eqn{\lambda_i / \lambda_1 \leq } \code{eig.tol}, default:
    \code{eig_tol = 0.01}.}
}
\details{
  
  The matrix \code{B} is assumed to be symmetric positive semi-definite 
  (negative eigenvalues not tolerated).  We start by computing the
  eigenvalue decomposition of \code{B}. Small eigenvalues (based on the
  \code{eig_tol} threshold) and their corresponding eigenvectors are
  then dropped. The result is the product of the remaining eigenvectors
  matrix (stored by columns) and the diagonal matrix formed by the
  square-root of the associated eigenvalues.

}
\value{
  
  If \code{B} is a \eqn{n \times n}{n * n}-matrix (symmetric
   positive semi-definite) and \eqn{\lambda_1,\dots,\lambda_k} its
  eigenvalues such as \eqn{\lambda_i / \lambda_1 > } \code{eig.tol}
  \eqn{,i=1,\dots,k} and let \code{UU} be the \eqn{n \times k}{n *
  k}-matrix formed by the associated eigenvectors ,i.e
  \eqn{UU=[U_1|\dots|U_k]}, the result factor matrix \code{C} is given
  by : \eqn{C=UU*diag(\sqrt{\lambda_i}), i=1,\dots,k}. Hence, we have
  \eqn{B=C*C^t}.

}
\author{
Tarak Kharrat   
}
\note{
  The result does not hold if the matrix is not symmetric or has too
  many non-positive eigenvalues. 
}
\seealso{
        \code{chol}
}
\examples{
sigma=matrix(c(1, 0.9, -0.3, 0.9, 1, -0.4, -0.3, -0.4, 1), ncol = 3)
## default eig_tol=0.01
mat=Robust_chol(sigma)
print(sum(mat \%*\%t(mat)-sigma))
## set eig_tol=0.1 to drop smallest eigenvalue
mat=Robust_chol(sigma,0.1)
print(sum(mat \%*\%t(mat)-sigma))
}
\keyword{Matrix-Factorization}

