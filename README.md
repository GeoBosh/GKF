This is a fast and flexible implementation of the Kalman filter, which can deal
with missing values and non-positive definite matrices for the variance of the
disturbances of the measurement equation. It is mostly written in C++ and relies
fully on linear algebra subroutines contained in the Armadillo library. Due to
the speed of the filter, the fitting of high-dimensional linear state space
models to large datasets is feasible. The package also treats nonlinear and
non-Gaussian models and allows signal smoothing and sampling from the (signal)
posterior distribution.


# Installing GKF

You can install the [development version](https://github.com/GeoBosh/GKF) of `GKF` from Github:

    library(devtools)
    install_github("GeoBosh/GKF")

See the [vignette](https://github.com/GeoBosh/GKF/blob/master/vignettes/vignette.pdf).

The package passes the CRAN quality control checks, if there is interest, we would publish it
there.

