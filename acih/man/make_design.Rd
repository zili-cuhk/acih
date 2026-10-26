\name{make_design}
\alias{make_design}
\alias{calibrate_design}
\alias{gen_data}
\alias{gen_X}
\alias{population_target}
\title{Simulation designs and data generation}
\description{
\code{make_design} defines a scenario (dimension, coefficients, covariate
design, censoring, baseline, horizon); \code{calibrate_design} calibrates the
censoring parameter on an independent Monte Carlo sample; \code{gen_data}
generates data from the additive hazards model with bounded covariates;
\code{gen_X} generates covariates; \code{population_target} computes the
least-false parameter on a fixed set from a large sample.
}
\value{A design list, a data list (\code{time}, \code{status}, \code{X}), a
covariate matrix, or a coefficient vector.}
