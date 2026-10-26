\name{acih-package}
\alias{acih-package}
\alias{acih}
\docType{package}
\title{Active-coordinate inference for sparse additive hazards regression under informative subsampling}
\description{
Implements the method and the numerical studies of the manuscript: Lin--Ying
estimating equations with risk-set centring, GSELO screening with EBIC on a
screening split, split-full and inverse-probability-weighted refits under
with-replacement and Poisson (including event-complete) subsampling with
inference-oriented allocation, and the total sandwich covariance.
}
\details{
Simulation study on a cluster: install the package, create a run directory
with \code{setup_run_dir("acih_run")}, then run \code{bash hpc/submit_all.sh}
there. See \code{\link{design_table}} for the designs and
\code{\link{ah_setup}} for the estimating equations.
}
