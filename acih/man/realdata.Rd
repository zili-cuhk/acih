\name{realdata}
\alias{realdata}
\alias{primary_split}
\alias{analyze_split}
\alias{ah_fit_survey}
\alias{km}
\alias{app_cohorts}
\alias{app_task_table}
\alias{load_cohort}
\alias{run_app_task}
\alias{summarize_apps}
\alias{fig_app_efficiency}
\alias{fig_app_km}
\title{Data-analysis protocol (Section 5)}
\description{
\code{primary_split} makes a 50/50 split; \code{analyze_split} screens one
split (GSELO with joint EBIC, or a fixed budget of the largest unpenalised
coefficients), refits the selected set on the inference split (split-full,
with-replacement IPW, Poisson designs and the case-cohort comparator at
several budgets) and times every step. \code{run_app_task} runs one split of
one cohort (an LSF array index; split 1 is the prespecified primary split,
which also carries the budget and survey sensitivity analyses);
\code{summarize_apps} collects all splits into CSV tables and figures.
\code{ah_fit_survey} is the fixed-support survey-weighted fit with stratified
PSU linearisation and \code{km} a base-R Kaplan-Meier estimator.
}
\usage{
primary_split(n, seed)
analyze_split(time, status, X, tau, split, screen = c("gselo", "topk"), budget = 10, cap = 20,
              r_vec = APP_R, seed_draw = 1, subsample = TRUE)
run_app_task(task, dir = ".", n_splits = 100, r_vec = NULL)
summarize_apps(dir = ".", cohorts = names(app_cohorts()))
}
\value{Lists of result tables; files are written under \code{dir/app}.}
