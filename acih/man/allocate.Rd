\name{allocate}
\alias{allocate}
\alias{ipw_probs}
\alias{draw_ipw}
\alias{poisson_inclusion}
\alias{poisson_scores}
\title{Subsampling probabilities}
\description{
\code{allocate} and \code{ipw_probs} give with-replacement probabilities for
the uniform, score-norm, Zuo et al. and inference-oriented rules, stratified by
event status and stabilized; \code{draw_ipw} draws indices by inverse CDF from
common uniforms. \code{poisson_scores} and \code{poisson_inclusion} give
Poisson inclusion probabilities with expected size \eqn{r}, including the
event-complete design.
}
\value{Probability vectors or sampled indices.}
