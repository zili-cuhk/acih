## Command-line entry point; run from a run directory created by acih::setup_run_dir().
##   simulation:   Rscript acih.R test | prepare | list | ntask | <task> | time [design ids] | summarize
##   applications: Rscript acih.R app_list | app_ntask | app <task> | app_summarize
##   plasmode:     Rscript acih.R plasmode_prepare | plasmode_ntask | plasmode <task> | plasmode_summarize
##   screening criterion check: Rscript acih.R screencmp_ntask | screencmp <task> | screencmp_summarize
##   figures only (from the summary CSV files): Rscript acih.R figures
suppressPackageStartupMessages(library(acih))
a <- commandArgs(trailingOnly = TRUE); if (!length(a)) a <- "list"
dir <- Sys.getenv("ACIH_DIR", getwd())
switch(a[1],
  test      = smoke_test(),
  prepare   = prepare_designs(dir),
  list      = { tt <- sim_task_table(load_designs(dir)); print(tt, row.names = FALSE); cat("tasks:", nrow(tt), "\n") },
  ntask     = cat(nrow(sim_task_table(load_designs(dir))), "\n", sep = ""),
  time      = time_designs(dir, if (length(a) > 1) a[-1] else NULL),
  summarize = summarize_results(dir),
  app_list  = { tt <- app_task_table(dir); print(table(tt$cohort)); cat("tasks:", nrow(tt), "\n") },
  app_ntask = cat(nrow(app_task_table(dir)), "\n", sep = ""),
  app       = run_app_task(as.integer(a[2]), dir),
  app_summarize = summarize_apps(dir),
  plasmode_prepare = prepare_plasmode(dir),
  plasmode_ntask = cat(nrow(plasmode_task_table(dir)), "\n", sep = ""),
  plasmode  = run_plasmode_task(as.integer(a[2]), dir),
  plasmode_summarize = summarize_plasmode(dir),
  screencmp_ntask = cat(nrow(screencmp_task_table(dir)), "\n", sep = ""),
  screencmp = run_screencmp_task(as.integer(a[2]), dir),
  screencmp_summarize = summarize_screencmp(dir),
  figures = make_figures_from_csv(dir),
  run_task(as.integer(a[1]), dir))
invisible(NULL)
