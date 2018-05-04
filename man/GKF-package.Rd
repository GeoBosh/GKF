\name{GKF-package}
\alias{GKF-package}
\alias{GKF}
\docType{package}
\title{
  \Sexpr[stage=build,results=hide]{pd <- packageDescription("GKF")}
  \Sexpr[stage=build,results=hide]{lb <- library(help="GKF", character.only=TRUE)}
  \Sexpr[stage=build]{pd$Title}
}
\description{
  \Sexpr[stage=build]{pd$Description}
}
\details{
   \tabular{ll}{
   Package:  \tab \Sexpr[stage=build]{pd$Package}\cr
   Type:     \tab \Sexpr[stage=build]{pd$Type}\cr
   Version:  \tab \Sexpr[stage=build]{pd$Version} \cr
   Date:     \tab \Sexpr[stage=build]{pd$Date}\cr
   License:  \tab \Sexpr[stage=build]{pd$License}\cr
   LazyLoad: \tab \Sexpr[stage=build]{pd$LazyLoad}\cr
   Built:    \tab \Sexpr[stage=build]{pd$Built}\cr
   }

   Index:
  \Sexpr[stage=build,results=rd]{paste("\\\\preformatted{", paste(if(!is.null(lb$info[[2]])) lb$info[[2]] else "No entries",collapse="\n"), "}", sep="")}

% ~~ An overview of how to use the package, including the most important ~~
% ~~ functions ~~
}
\author{
  \Sexpr[stage=build]{pd$Author}

Maintainer: \Sexpr[stage=build]{pd$Maintainer}
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
\keyword{ package }
\keyword{Kalman-Filter }
\seealso{
  %~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
   \code{dlm}, \code{KFAS}, \code{FKF} 
}
\examples{
"See Examples in package-vignette"
}
