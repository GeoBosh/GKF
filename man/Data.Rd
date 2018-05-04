\name{GKF_data}
\alias{GKF_data}
\alias{GKF.data}
\docType{data}
\title{
  Data used in the vignette examples
}
\description{
  The package is shipped with 4 data sets used in the examples described in the
  vignette.
}
\usage{data(GKF_data)}
\format{
  The format is:
  \preformatted{  
List of 4
 $ longGoals : num [1:10, 1:180] 0 0 1 0 1 0 1 2 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:10] "Arsenal" "Aston Villa" "Chelsea" "Everton" ...
  .. ..$ : NULL
 $ shortGoals: num [1:10, 1:36] 0 0 1 2 0 2 5 3 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:10] "Arsenal" "Aston Villa" "Chelsea" "Everton" ...
  .. ..$ : NULL
 $ invest    :'data.frame':	41 obs. of  2 variables:
  ..$ V1: num [1:41] 44.7 45.9 50.2 44 55.7 ...
  ..$ V2: num [1:41] 2617 3417 4173 4566 4822 ...
 $ CAPMData  : mts [1:120, 1:6] -4.6 -1.7 4.9 7.7 -1.1 -4.3 2.8 5.6 6.4 -6.9 ...
  ..- attr(*, "tsp")= num [1:3] 1978 1988 12
  ..- attr(*, "class")= chr [1:3] "mts" "ts" "matrix"
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:6] "MOBIL" "IBM" "WEYER" "CITCRP" ...
  }
}
\details{
  
  The first two data sets \code{longGoals} and \code{shortGoals} are
  goals scored by 10 English Premier league in the past seasons. Only
  the length of the observation time is different.  Those data sets are
  used in the nonlinear non-Gaussian example.

  The \code{invest} corresponds to the Seemingly unrelated time series
  example proposed in Petris (2009, chap 3.3.2) and representing the
  annual Denmark and Spain investments .

  The last data set \code{CAPMData} corresponds to the Seemingly
  unrelated regression models example proposed in Petris (2009, chap
  3.3.3) and representing the the monthly returns for the Mobil, IBM,
  Weyer,and Citicorp stocks, 1978.1-1987.12.

}
\source{

  The first two data sets were collected and adapted from
  \emph{http://www.football-data.co.uk/}.  The last two were downloaded
  from Petris (2009) book's website
  \emph{http://definetti.uark.edu/~gpetris/dlm/}.

}
\references{
        % bibentry: petris2009dynamic
Petris G, Petrone S and Campagnoli P (2009).
\emph{Dynamic linear models with R}.
Springer.
% end:bibentry:  petris2009dynamic
}
\examples{
data(GKF_data)
str(GKF.data)
}
\keyword{datasets}
