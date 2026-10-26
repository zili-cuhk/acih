\name{design_table}
\alias{design_table}
\alias{prepare_designs}
\alias{load_designs}
\alias{sim_task_table}
\alias{run_task}
\alias{time_designs}
\alias{summarize_results}
\alias{setup_run_dir}
\alias{smoke_test}
\alias{ipw_methods_default}
\alias{run_rep_known}
\alias{run_rep_pipeline}
\alias{run_rep_marginal}
\title{Simulation study of the manuscript}
\description{
\code{design_table} lists the designs of Studies 1--3; \code{prepare_designs}
calibrates them once into \code{dir/out/designs.rds}; \code{sim_task_table}
splits them into tasks; \code{run_task} runs one task into \code{dir/out/};
\code{time_designs} times one replicate per design; \code{summarize_results}
writes the tables and figures into \code{dir/results/}; \code{setup_run_dir}
copies the command-line script, LSF scripts and application scripts into a
run directory; \code{smoke_test} runs each study type once at reduced size.
\code{run_rep_known}, \code{run_rep_pipeline} and \code{run_rep_marginal} run
one replicate of the refit-stage, complete-pipeline and marginal-route
studies, with the methods of \code{ipw_methods_default}.
}
\usage{
design_table()
prepare_designs(dir = ".", n_cal = 2e5, n_big = 1e6, verbose = TRUE)
load_designs(dir = ".")
sim_task_table(designs)
run_task(task, dir = ".")
time_designs(dir = ".", ids = NULL)
summarize_results(dir = ".")
setup_run_dir(dir)
smoke_test()
}
\arguments{
\item{dir}{run directory.}
\item{n_cal}{Monte Carlo size for censoring calibration.}
\item{n_big}{Monte Carlo size for least-false targets.}
\item{verbose}{print calibration diagnostics.}
\item{designs}{output of \code{load_designs}.}
\item{task}{task index.}
\item{ids}{design identifiers.}
}
\value{Mostly called for side effects.}
