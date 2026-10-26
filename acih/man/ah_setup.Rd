\name{ah_setup}
\alias{ah_setup}
\alias{ah_bV}
\alias{q_scores}
\alias{xi_scores}
\alias{ah_fit_full}
\alias{ah_fit_ipw}
\alias{ah_fit_poisson}
\alias{ah_fit_casecohort}
\alias{pilot_fit}
\title{Lin--Ying estimating equations, refits and sandwich covariances}
\description{
\code{ah_bV} computes the ambient-dimension \eqn{b} and \eqn{V} blockwise for
screening. \code{ah_setup} computes the per-subject quantities on a selected
set (\eqn{e_i}, \eqn{G_i}, Breslow residual parts) on an inference split;
\code{q_scores} and \code{xi_scores} return the score contributions
\eqn{q_i(\beta)} and the martingale residuals \eqn{\xi_i(\beta)}.
\code{ah_fit_full} is the split-full refit with the martingale sandwich;
\code{ah_fit_ipw} and \code{ah_fit_poisson} are the with-replacement and
Poisson IPW refits with the total sandwich
\eqn{H^{-1}(\Omega/n+\Gamma/r)H^{-T}}; \code{pilot_fit} fits the pilot on the
screening split. \code{ah_fit_casecohort} is the case-cohort estimator of
Kulich and Lin (2000, Biometrika 87, 73--87) with their variance estimator,
used as a comparator on the event-complete subsample.
}
\usage{
ah_setup(time, status, XA, tau, w = NULL)
ah_bV(time, status, X, tau, block = 5000L)
}
\arguments{
\item{time, status}{observed times and event indicators.}
\item{XA, X}{covariates of the selected set, or all covariates.}
\item{tau}{end of the analysis window.}
\item{w}{optional analysis weights.}
\item{block}{row block size.}
}
\value{Lists of the quantities described above; see the source for details.}
