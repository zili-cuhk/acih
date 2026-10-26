\name{plasmode}
\alias{plasmode}
\alias{plasmode_design}
\alias{plasmode_data}
\alias{run_rep_plasmode}
\alias{run_rep_plasmode_known}
\alias{prepare_plasmode}
\alias{plasmode_task_table}
\alias{run_plasmode_task}
\alias{summarize_plasmode}
\title{Plasmode simulation built on a real cohort}
\description{
\code{plasmode_design} keeps the covariates of a processed cohort, fits the
reference support on the whole cohort, and sets up an additive hazards model
with a positive constant baseline and uniform censoring calibrated to the
observed event proportion within the horizon; \code{plasmode_data} draws new
outcomes. \code{run_rep_plasmode_known} refits a fixed reference set
containing the true support on a new 50/50 split (split-full, with-replacement
and Poisson designs, the case-cohort comparator, and the inference-split-only
and subsampling-only variance comparators); \code{run_rep_plasmode} runs the
complete pipeline including GSELO screening. \code{prepare_plasmode},
\code{plasmode_task_table}, \code{run_plasmode_task} and
\code{summarize_plasmode} drive the study from the results of the data
analysis (true support: variables selected in at least 95 of 100 splits;
reference set: the primary-split selected set).
}
\usage{
plasmode_design(dat, covs, support, tau, floor_rate = NULL, n_cal = 2e5, seed = 1)
plasmode_data(pd)
run_rep_plasmode_known(pd, A, seed, r_vec,
  methods = ipw_methods_default()[c("UIPW", "SIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW")])
prepare_plasmode(dir = ".", cohorts = names(app_cohorts()), stable = 0.95)
run_plasmode_task(task, dir = ".")
summarize_plasmode(dir = ".")
}
\value{Lists of replicate results and summary tables written under \code{dir/app}.}
