\name{InvMat}
\alias{InvMat}

\title{
        Smoothed Matrix Inverse
}
\description{

  The function computes the inverse of nearly singular matrices by
  inflating the diagonal (when needed).

}
\usage{
InvMat(input)
}
\arguments{
  \item{input}{A \code{matrix} to be inverted.}
}
\details{

  The matrix inverse is computed directly by the \code{inv} routine in
  \code{Armadillo} unless the condition number is smaller than
  \code{1e-12}. In this case, the matrix diagonal is inflated by small
  number (\code{1e-6}) and \code{inv} is used again. This trick avoids
  the use of the Cholesky decomposition which requires the matrix to be
  definite positive (and fails if the matrix is nearly singular).

}
\author{
  Tarak Kharrat
}
\seealso{
\code{solve}
}
\examples{
set.seed(345) 
M <- matrix(rnorm(16), 4, 4)
sum(solve(M) - InvMat(M))
}
\keyword{algebra}
\keyword{array}
