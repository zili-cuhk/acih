## 09_run.R -- task table and task runner for the simulation study, run directory set-up.
## All files live under a run directory `dir`:
##   dir/out/designs.rds, dir/out/<study>/<design>/reps_<first>_<last>.rds, dir/results/.

## replicates per task: chosen so that a task stays well below the 24-hour limit
sim_chunk <- function(des) {
  if (des$run %in% c("known", "marginal")) return(50)
  if (des$n >= 200000 || des$p >= 1000) return(2)
  5
}
BASE_SEED <- 20260711L

load_designs <- function(dir = ".") {
  f <- file.path(dir, "out", "designs.rds")
  if (!file.exists(f)) stop("no designs found in ", f, "; run prepare_designs() first")
  readRDS(f)
}

sim_task_table <- function(designs) {
  tt <- do.call(rbind, lapply(seq_along(designs), function(i) {
    des <- designs[[i]]; ch <- sim_chunk(des); first <- seq(1, des$R, by = ch)
    data.frame(task = NA_integer_, design_index = i, id = names(designs)[i], study = des$study,
               first = first, last = pmin(first + ch - 1, des$R), stringsAsFactors = FALSE)
  }))
  tt$task <- seq_len(nrow(tt)); tt
}

run_one <- function(des, seed) {
  switch(des$run,
    known = run_rep_known(des, des$extra$m, des$extra$r, seed, target = des$lf_target),
    pipeline = run_rep_pipeline(des, seed, r_vec = des$extra$r %||% c(1000, 5000),
                                methods = ipw_methods_default()[c("UIPW", "SIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW")],
                                naive = isTRUE(des$extra$naive)),
    marginal = run_rep_marginal(des, seed))
}

## Run one task (an LSF array index).  Existing output is not recomputed.
run_task <- function(task, dir = ".") {
  designs <- load_designs(dir); tt <- sim_task_table(designs)
  if (task < 1 || task > nrow(tt)) stop("no such task: ", task)
  tk <- tt[task, ]; des <- designs[[tk$id]]
  od <- file.path(dir, "out", tk$study, tk$id); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  of <- file.path(od, sprintf("reps_%04d_%04d.rds", tk$first, tk$last))
  if (file.exists(of)) { message("exists: ", of); return(invisible(of)) }
  res <- list()
  for (r in tk$first:tk$last) {
    t0 <- tic()
    res[[as.character(r)]] <- tryCatch(run_one(des, BASE_SEED + 10000L * tk$design_index + r),
                                       error = function(e) list(error = conditionMessage(e)))
    message(sprintf("[%s] rep %d %.1fs", tk$id, r, tic() - t0))
  }
  saveRDS(list(id = tk$id, reps = res, R_version = R.version.string), of)
  invisible(of)
}

## Time one replicate of selected designs (to check the chunk sizes against the 24-hour limit).
time_designs <- function(dir = ".", ids = NULL) {
  designs <- load_designs(dir)
  if (is.null(ids)) ids <- names(designs)
  out <- data.frame(id = ids, seconds = NA_real_, chunk = NA_real_, task_hours = NA_real_)
  for (k in seq_along(ids)) {
    des <- designs[[ids[k]]]; t0 <- tic()
    invisible(run_one(des, BASE_SEED + 999L))
    out$seconds[k] <- tic() - t0; out$chunk[k] <- sim_chunk(des)
    out$task_hours[k] <- out$seconds[k] * out$chunk[k] / 3600
    message(sprintf("%-22s %8.1f s per replicate; %d per task = %.2f h", ids[k], out$seconds[k],
                    out$chunk[k], out$task_hours[k]))
  }
  invisible(out)
}

## Copy the command-line script, LSF scripts and application scripts into a run directory.
setup_run_dir <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (s in c("logs", "out", "results", "hpc", "app")) dir.create(file.path(dir, s), showWarnings = FALSE)
  file.copy(system.file("scripts", "acih.R", package = "acih"), dir, overwrite = TRUE)
  file.copy(list.files(system.file("hpc", package = "acih"), full.names = TRUE), file.path(dir, "hpc"), overwrite = TRUE)
  file.copy(list.files(system.file("app", package = "acih"), full.names = TRUE), file.path(dir, "app"), overwrite = TRUE)
  for (s in c("data/msk_met_2021", "data/gene_panels", "data/nhis", "out", "results"))
    dir.create(file.path(dir, "app", s), showWarnings = FALSE, recursive = TRUE)
  message("run directory ready: ", normalizePath(dir))
  invisible(dir)
}

## One replicate of each study type at reduced size.
smoke_test <- function() {
  t0 <- tic()
  dk <- calibrate_design(make_design(n = 1, p = 8, cens = "uniform", cens_rate = 0.5), n_cal = 2e4)
  r1 <- run_rep_known(dk, 3000, c(300, 1000), seed = 1)
  stopifnot(!r1$gate_fail, nrow(r1$res) > 0)
  dp <- calibrate_design(make_design(n = 6000, p = 60, cens = "uniform", cens_rate = 0.5), n_cal = 2e4)
  r2 <- run_rep_pipeline(dp, seed = 2, r_vec = 500, naive = TRUE,
                         methods = ipw_methods_default()[c("UIPW", "IIPW", "PUIPW", "PIIPW")])
  stopifnot(!r2$gate_fail, nrow(r2$res) > 0)
  dm <- calibrate_design(make_design(n = 600, p = 300, s_design = "reinf4", cens = "uniform", cens_rate = 0.25), n_cal = 2e4)
  r3 <- run_rep_marginal(dm, seed = 3)
  stopifnot(nrow(r3$res) == 2)
  message(sprintf("smoke test passed in %.1f s", tic() - t0))
  invisible(TRUE)
}


## ---- screening-criterion comparison (EBIC with log n versus log events) ----
## Screening only: same data, split and GSELO path as the pipeline replicates
## (same seeds), with the BIC term |A| log n_scr or |A| log N_event.
SCREENCMP_R <- 200L
screencmp_designs <- function(dir = ".") {
  designs <- load_designs(dir)
  designs[vapply(designs, function(d) identical(d$run, "pipeline"), TRUE)]
}
screencmp_task_table <- function(dir = ".") {
  designs <- load_designs(dir); ids <- names(screencmp_designs(dir))
  tt <- do.call(rbind, lapply(ids, function(id) {
    des <- designs[[id]]; R <- min(des$R, SCREENCMP_R)
    ch <- if (des$n >= 200000 || des$p >= 1000) 5L else 10L; first <- seq(1, R, by = ch)
    data.frame(design_index = match(id, names(designs)), id = id, first = first,
               last = pmin(first + ch - 1, R), stringsAsFactors = FALSE)
  }))
  tt$task <- seq_len(nrow(tt)); tt
}
screencmp_rep <- function(des, seed) {
  set.seed(seed)
  d <- gen_data(des); n <- length(d$time)
  if (mean(d$time >= des$tau) < 0.05) return(NULL)
  scr_idx <- sort(sample.int(n, floor(n / 2)))
  tm <- d$time[scr_idx]; stt <- d$status[scr_idx]; X <- d$X[scr_idx, , drop = FALSE]
  fac <- time_rate(tm, stt, des$tau); tm <- tm * fac; tau <- des$tau * fac
  bv <- ah_bV(tm, stt, X, tau); lam <- log(max(bv$n_event, 1)) / bv$n
  gp <- gselo_path(bv$V, bv$b, lam, theta_grid_default())
  S0 <- des$S0
  do.call(rbind, lapply(c(n = log(bv$n), events = log(max(bv$n_event, 2))), function(lp) {
    sel <- ebic_select(bv$V, bv$b, gp, bv$n, 0.25, eta_grid_default(), 20, ncol(X), loss_mult = 2, log_pen = lp)
    A <- sel$A
    data.frame(size = length(A), sure = all(S0 %in% A), exact = setequal(A, S0),
               n_fp = length(setdiff(A, S0)), n_event = bv$n_event)
  })) -> o
  o$criterion <- c("log n", "log events"); o
}
run_screencmp_task <- function(task, dir = ".") {
  designs <- load_designs(dir); tt <- screencmp_task_table(dir); tk <- tt[task, ]; des <- designs[[tk$id]]
  od <- file.path(dir, "out", "screencmp", tk$id); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  of <- file.path(od, sprintf("reps_%04d_%04d.rds", tk$first, tk$last))
  if (file.exists(of)) return(invisible(of))
  res <- lapply(tk$first:tk$last, function(r) {
    o <- tryCatch(screencmp_rep(des, BASE_SEED + 10000L * tk$design_index + r), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o })
  saveRDS(do.call(rbind, res), of); invisible(of)
}
summarize_screencmp <- function(dir = ".") {
  ids <- list.dirs(file.path(dir, "out", "screencmp"), recursive = FALSE)
  out <- do.call(rbind, lapply(ids, function(d) {
    x <- do.call(rbind, lapply(list.files(d, full.names = TRUE), readRDS))
    if (is.null(x) || !nrow(x)) return(NULL)
    a <- stats::aggregate(cbind(sure, exact, size, n_fp, n_event) ~ criterion, x, mean)
    a$reps <- as.vector(table(x$criterion)[a$criterion]); a$design <- basename(d); a
  }))
  dir.create(file.path(dir, "results"), showWarnings = FALSE)
  utils::write.csv(out, file.path(dir, "results", "screencmp.csv"), row.names = FALSE)
  print(out, row.names = FALSE, digits = 3); invisible(out)
}
